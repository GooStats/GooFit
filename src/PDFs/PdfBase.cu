#include "goofit/BinnedDataSet.h"
#include "goofit/FitControl.h"
#include "goofit/PdfBase.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/Variable.h"

// This is code that belongs to the PdfBase class, that is,
// it is common across all implementations. But it calls on device-side
// functions, and due to the nvcc translation-unit limitations, it cannot
// sit in its own object file; it must go in the CUDAglob.cu. So it's
// off on its own in this inline-cuda file, which GooPdf.cu
// should include.

__host__ void PdfBase::copyParams(const std::vector<double> &pars) {
  // copyParams method performs eponymous action!

  for (unsigned int i = 0; i < pars.size(); ++i) {
    host_params[i] = pars[i];

    if (std::isnan(host_params[i])) {
      std::cout << " agh, parameter is NaN, die " << i << std::endl;
      abortWithCudaPrintFlush(__FILE__, __LINE__, "NaN in parameter");
    }
  }

  MEMCPY_TO_SYMBOL(cuda_array, host_params, pars.size() * sizeof(fptype), 0, cudaMemcpyHostToDevice);
}

__host__ void PdfBase::copyParams() const {
  // Copies values of Variable objects
  parCont pars;
  getParameters(pars);
  int maxIndex = -1;
  for (auto &par: pars) {
    if (maxIndex < par->getIndex()) {
      maxIndex = par->getIndex();
    }
  }
  std::vector<double> values(host_params, host_params + maxIndex + 1);
  for (auto &par: pars) {
    values[par->getIndex()] = par->value;
  }
  copyParams(values);
}

__host__ void PdfBase::copyNormFactors() {
  //MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice);
  SYNCH();// Ensure normalisation integrals are finished
}

__host__ void PdfBase::initialiseIndices(std::vector<unsigned int> pindices) {
  // Structure of the individual index array: Number of parameters, then the indices
  // requested by the subclass (which will be interpreted by the subclass kernel),
  // then the number of observables, then the observable indices. Notice that the
  // observable indices are not set until 'setIndices' is called, usually from setData;
  // here we only reserve space for them by setting totalParams.
  // This is to allow index sharing between PDFs - all the PDFs must be constructed
  // before we know what observables exist.

  if (totalParams + pindices.size() >= maxIndicies) {
    std::cout << "Major problem with pindices size: " << totalParams << " + " << pindices.size() << " >= " << maxIndicies << std::endl;
  }

  assert(totalParams + pindices.size() < maxIndicies);
  host_indices[totalParams] = pindices.size();
  for (unsigned int i = 1; i <= host_indices[totalParams]; ++i) {
    host_indices[totalParams + i] = pindices[i - 1];
  }
  host_indices[totalParams + pindices.size() + 1] = observables.size();

  parameters = totalParams;
  totalParams += (2 + pindices.size() + observables.size());
  pdfName[parameters] = getName();

  /*
  std::cout << " | "
	    << parameters << " "
	    << totalParams << " "
	    << cuda_array << " "
	    << paramIndices << " "
	    << std::endl;
  */
  MEMCPY_TO_SYMBOL(paramIndices, host_indices, totalParams * sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
}

__host__ void PdfBase::setData(std::vector<std::map<Variable *, fptype>> &__attribute__((__unused__)) data) {
  abortWithCudaPrintFlush(__FILE__, __LINE__, "not possible to use this method. one Variable could correspond to different index, to save space");
}

__host__ void PdfBase::recursiveSetIndices() {
  for (auto &component: components) {
    component->recursiveSetIndices();
  }

  int numParams = host_indices[parameters];
  int counter = 0;
  for (auto v = obsBegin(); v != obsEnd(); ++v) {
    host_indices[parameters + 2 + numParams + counter] = counter;
    counter++;
  }
  generateNormRange();
}

__host__ void PdfBase::setIndices() {
  recursiveSetIndices();
  MEMCPY_TO_SYMBOL(paramIndices, host_indices, totalParams * sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
}

__host__ void PdfBase::setData(UnbinnedDataSet *data) {
  if (pdfId == -1) {
    registerPdf();
  }
  if (dev_event_array[pdfId] != nullptr) {
    gooFree(dev_event_array[pdfId]);
    SYNCH();
    dev_event_array[pdfId] = nullptr;
  }

  setIndices();
  int dimensions = observables.size();
  numEntries = data->getNumEvents();
  numEvents = numEntries;
  if (fitControl->binnedFit()) {
    setFitControl(new UnbinnedNllFit());
  }

  if (numEntries > 0) {
    auto *host_array = new fptype[numEntries * dimensions];
    for (unsigned int i = 0; i < numEntries; ++i) {
      int j = 0;
      for (auto v = obsBegin(); v != obsEnd(); ++v) {
        fptype currVal = data->getValue((*v), i);
        int position = host_indices[parameters + 2 + host_indices[parameters] + j];
        host_array[i * dimensions + position] = currVal;
        ++j;
      }
    }

    gooMalloc(reinterpret_cast<void **>(&(dev_event_array[pdfId])), dimensions * numEntries * sizeof(fptype));
    MEMCPY(dev_event_array[pdfId], host_array, dimensions * numEntries * sizeof(fptype), cudaMemcpyHostToDevice);
    SYNCH();
    MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice);

    delete[] host_array;
  }
}

__host__ void PdfBase::setData(BinnedDataSet *data) {
  if (pdfId == -1) {
    registerPdf();
  }
  if (dev_event_array[pdfId] != nullptr) {
    gooFree(dev_event_array[pdfId]);
    dev_event_array[pdfId] = nullptr;
  }

  setIndices();
  numEvents = 0;
  numEntries = data->getNumBins();
  if (numEntries <= 0) {
    abortWithCudaPrintFlush(__FILE__, __LINE__, "0 entries. check the numbins of the variable of your data set", this);
  }
  int dimensions = 2 + observables.size();// Bin center (x,y, ...), bin value, and bin volume.
  if (!fitControl->binnedFit()) {
    setFitControl(new BinnedNllFit());
  }

  if (numEntries > 0) {
    auto *host_array = new fptype[numEntries * dimensions];

    for (unsigned int i = 0; i < numEntries; ++i) {
      int j = 0;
      for (auto v = obsBegin(); v != obsEnd(); ++v) {
        int position = host_indices[parameters + 2 + host_indices[parameters] + j];
        assert(position == j);
        host_array[i * dimensions + position] = data->getBinCenter((*v), i);
        ++j;
      }

      host_array[i * dimensions + observables.size() + 0] = data->getBinContent(i);
      host_array[i * dimensions + observables.size() + 1] = fitControl->binErrors() ? data->getBinError(i) : data->getBinVolume(i);
      numEvents += data->getBinContent(i);
    }

    gooMalloc(reinterpret_cast<void **>(&(dev_event_array[pdfId])), dimensions * numEntries * sizeof(fptype));
    MEMCPY(dev_event_array[pdfId], host_array, dimensions * numEntries * sizeof(fptype), cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(functorConstants, &numEvents, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    delete[] host_array;
  }
}

__host__ void PdfBase::generateNormRange() {
  if (normRanges != nullptr) {
    gooFree(normRanges);
  }
  gooMalloc(reinterpret_cast<void **>(&normRanges), 3 * observables.size() * sizeof(fptype));

  auto *host_norms = new fptype[3 * observables.size()];
  int counter = 0;// Don't use index in this case to allow for, eg,
  // a single observable whose index is 1; or two observables with indices
  // 0 and 2. Make one array per functor, as opposed to variable, to make
  // it easy to pass MetricTaker a range without worrying about which parts
  // to use.
  for (auto v = obsBegin(); v != obsEnd(); ++v) {
    host_norms[3 * counter + 0] = (*v)->lowerlimit;
    host_norms[3 * counter + 1] = (*v)->upperlimit;
    host_norms[3 * counter + 2] = integrationBins > 0 ? integrationBins : (*v)->numbins;
    counter++;
  }

  MEMCPY(normRanges, host_norms, 3 * observables.size() * sizeof(fptype), cudaMemcpyHostToDevice);
  delete[] host_norms;
}

void PdfBase::clearCurrentFit() const {
  totalParams = 0;
  gooFree(dev_event_array[pdfId]);
  dev_event_array[pdfId] = nullptr;
}

__host__ void PdfBase::printProfileInfo(bool __attribute__((__unused__)) topLevel) {
#ifdef PROFILING
  if (topLevel) {
    cudaError_t err = MEMCPY_FROM_SYMBOL(host_timeHist, timeHistogram, 10000 * sizeof(fptype), 0, cudaMemcpyDeviceToHost);
    if (cudaSuccess != err) {
      std::cout << "Error on copying timeHistogram: " << cudaGetErrorString(err) << std::endl;
      return;
    }

    std::cout << getName() << " : " << getFunctionIndex() << " " << host_timeHist[100 * getFunctionIndex() + getParameterIndex()] << std::endl;
    for (unsigned int i = 0; i < components.size(); ++i) {
      components[i]->printProfileInfo(false);
    }
  }
#endif
}


gooError gooMalloc(void **target, size_t bytes) {
  // Thrust 1.7 will make the use of THRUST_DEVICE_BACKEND an error
#if THRUST_DEVICE_SYSTEM != THRUST_DEVICE_SYSTEM_CUDA
  target[0] = malloc(bytes);
  if (target[0] != nullptr) {
    return gooSuccess;
  }
  return gooErrorMemoryAllocation;
#else
  if (cudaMalloc(target, bytes) != cudaSuccess)
    throw std::runtime_error("cannot allocate enough memory on the device");
  return (gooError) cudaSuccess;
#endif
}

gooError gooFree(void *ptr) {
  gooError ret;
  // Thrust 1.7 will make the use of THRUST_DEVICE_BACKEND an error
#if THRUST_DEVICE_SYSTEM != THRUST_DEVICE_SYSTEM_CUDA
  free(ptr);
  ret = gooSuccess;
#else
  ret = (gooError) cudaFree(ptr);
#endif
  ptr = nullptr;
  return ret;
}
