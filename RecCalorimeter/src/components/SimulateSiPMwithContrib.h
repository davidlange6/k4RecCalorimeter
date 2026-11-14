#ifndef SimulateSiPMwithContrib_h
#define SimulateSiPMwithContrib_h 1

#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/TimeSeriesCollection.h"

#include "k4FWCore/Transformer.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "TRandom3.h"

#include "GaudiKernel/ToolHandle.h"

// Check for SiPMSensor header location (similar to how DigiSiPM handles this)
#if __has_include("SiPMSensor.h")
#include "SiPMSensor.h"
#elif __has_include("sipm/SiPMSensor.h")
#include "sipm/SiPMSensor.h"
#endif

/** @class SimulateSiPMwithContrib
 *
 *  This algorithm was created to apply the SiPM response over the
 *  hit created with IDEA_o2 simulation.
 *  Contributions from such hits containt the photo-electrons (fired cells)
 *  and their time of arrivals at SiPMs.
 *  This algorithm takes all the input collections form the IDEA_o2
 *  hadronic calorimeter and returns corresponding digitized hits.
 *
 *  The algorithm simulates SiPM response including effects like:
 *  - Dark count rate (DCR)
 *  - Crosstalk (Xt)
 *  - Afterpulsing
 *  - Cell recovery time
 *  - Signal rise and fall times
 *
 *  @author Lorenzo Pezzotti
 *  @date   2026-02-17
 */

using SimulateSiPMwithContribOutput =
    std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>;

class SimulateSiPMwithContrib final
    : public k4FWCore::MultiTransformer<SimulateSiPMwithContribOutput(const edm4hep::SimCalorimeterHitCollection&,
                                                                      const edm4hep::EventHeaderCollection&)> {
public:
  SimulateSiPMwithContrib(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  SimulateSiPMwithContribOutput operator()(const edm4hep::SimCalorimeterHitCollection& simHits,
                                           const edm4hep::EventHeaderCollection& headers) const override;

private:
  // Per-event engine for the draws made by this algorithm (photon statistics, wavelength,
  // absorption, scintillation time). Reseeded from UniqueIDGenSvc at every event so that the
  // result for an event does not depend on the other events processed in the job.
  inline static thread_local TRandom3 m_engine;

  // readout name and segmentation (of specific type)
  // Not used for the moment, might be used in future
  // Gaudi::Property<std::string> m_readoutName{this, "readoutName", "", "name of the readout"};

  // unique ID service, used to seed the SiPM simulation for every event
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  std::unique_ptr<sipm::SiPMSensor> m_sensor;

  // Hamamatsu S14160-1310PS
  // Signal properties
  Gaudi::Property<double> m_sigLength{this, "signalLength", 200., "signal length in ns"};
  Gaudi::Property<double> m_sampling{this, "sampling", 0.1, "SiPM sampling rate in ns"};
  Gaudi::Property<double> m_risetime{this, "risetime", 1., "signal rise time in ns"};
  Gaudi::Property<double> m_falltimeFast{this, "falltimeFast", 6.5, "signal fast component decay time in ns"};

  // SiPM physical properties
  Gaudi::Property<double> m_sipmSize{this, "SiPMsize", 1.3, "width of photosensitive area in mm"};
  Gaudi::Property<double> m_cellPitch{this, "cellpitch", 10., "SiPM cell size in um"};
  Gaudi::Property<double> m_recovery{this, "recovery", 10.,
                                     "SiPM cell recovery time in ns"}; // https://arxiv.org/abs/2001.10322

  // Noise parameters
  Gaudi::Property<double> m_Dcr{this, "DCR", 120e3, "SiPM DCR"}; // dark count rate
  Gaudi::Property<double> m_Xt{this, "Xtalk", 0.01, "SiPM crosstalk"};
  Gaudi::Property<double> m_afterpulse{this, "afterpulse", 0.03, "afterpulse probability"};
  Gaudi::Property<double> m_snr{this, "SNR", 20., "SNR value in dB"}; // signal-to-noise ratio

  // integration parameters
  Gaudi::Property<double> m_gateStart{this, "gateStart", 5., "Integration gate starting time in ns"};
  Gaudi::Property<double> m_gateL{this, "gateLength", 95.,
                                  "Integration gate length in ns"}; // Should be approx 5 times fallTimeFast (see above)
  Gaudi::Property<double> m_thres{this, "threshold", 1.5,
                                  "Integration threshold"}; // Threshold in p.e. (1.5 to suppress DCR)

  // other parameters (attention, will override above parameters if set)
  Gaudi::Property<std::map<std::string, double>> m_params{this, "params", {}, "optional parameters"};

  // switch to select if you are processing scintillation or Cherenkov hits
  Gaudi::Property<bool> m_isCherenkov{this, "isCherenkov", false,
                                      "switch to select if you are processing scintillation or Cherenkov hits"};
  // SiPM efficiency
  Gaudi::Property<std::vector<double>> m_wavelen{
      this, "wavelength", {1000., 100.}, "wavelength vector in nm (decreasing order)"};
  Gaudi::Property<std::vector<double>> m_sipmEff{this, "sipmEfficiency", {1.0, 1.0}, "SiPM efficiency vs wavelength"};

  Gaudi::Property<double> m_scintDecaytime{this, "scintDecaytime", 2.8, "scintillation decay time in ns"};

  // option to store full waveform (for debugging)
  // Not used for the moment, might be used in future for debugging
  // Gaudi::Property<bool> m_storeFullWaveform{this, "storeFullWaveform", false, "Store full waveform for debugging"};

  // Suppression threshold parameters
  Gaudi::Property<double> m_suppressionThreshold{this, "suppressionTreshold", 0.001,
                                                 "Threshold on ADC integral to reject digi hit"};
};

#endif
