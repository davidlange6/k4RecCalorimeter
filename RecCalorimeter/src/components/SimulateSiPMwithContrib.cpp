#include "SimulateSiPMwithContrib.h"
#include "DD4hep/DD4hepUnits.h"
#include <cmath>

DECLARE_COMPONENT(SimulateSiPMwithContrib)

SimulateSiPMwithContrib::SimulateSiPMwithContrib(const std::string& aName, ISvcLocator* aSvcLoc)
    : MultiTransformer(aName, aSvcLoc,
                       {KeyValue("inputHitCollection", "SimCalorimeterHitCollection"),
                        KeyValue("inputHeaderCollection", "EventHeader")},
                       {KeyValue("outputHitCollection", "CalorimeterHitCollection"),
                        KeyValue("outputLinkCollection", "CaloHitSimCaloHitLinkCollection")}) {}

StatusCode SimulateSiPMwithContrib::initialize() {
  StatusCode sc = MultiTransformer::initialize();

  if (sc.isFailure())
    return sc;

  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  // initialize SiPM properties
  sipm::SiPMProperties properties;
  properties.setSignalLength(m_sigLength);
  properties.setSize(m_sipmSize);
  properties.setDcr(m_Dcr);
  properties.setXt(m_Xt);
  properties.setSampling(m_sampling);
  properties.setRecoveryTime(m_recovery);
  properties.setPitch(m_cellPitch);
  properties.setAp(m_afterpulse);
  properties.setFallTimeFast(m_falltimeFast);
  properties.setRiseTime(m_risetime);
  properties.setSnr(m_snr);
  properties.setPdeType(sipm::SiPMProperties::PdeType::kSimplePde); // this sets PDE to a single number
  properties.setPde(1.0);
  // properties.setPdeSpectrum(m_wavelen, m_sipmEff); // if you want wavelen dependency

  // set other parameters if provided
  for (const auto& [key, value] : m_params.value())
    properties.setProperty(key, value);

  m_sensor = std::make_unique<sipm::SiPMSensor>(properties); // must be constructed from SiPMProperties

  info() << "SimulateSiPMwithEdep initialized" << endmsg;
  info() << properties << endmsg; // sipm::SiPMProperties has std::ostream& operator<<

  return StatusCode::SUCCESS;
}

SimulateSiPMwithContribOutput SimulateSiPMwithContrib::operator()(const edm4hep::SimCalorimeterHitCollection& simHits,
                                                                  const edm4hep::EventHeaderCollection& headers) const {
  // Seed the SimSiPM generator from (run, event, algorithm name), so that the result does not depend
  // on the processing order of the events and is reproducible from job to job
  const uint64_t seed = m_uidSvc->getUniqueID(headers, name());
  m_sensor->rng().seed(seed != 0 ? seed : 1); // a zero seed is not allowed by SimSiPM
  const uint32_t seed32 = static_cast<uint32_t>(seed ^ (seed >> 32));
  m_engine.SetSeed(seed32 != 0 ? seed32 : 1); // TRandom3 treats 0 as "seed from the clock"

  edm4hep::CalorimeterHitCollection digiHits;
  edm4hep::CaloHitSimCaloHitLinkCollection links;

  // loop through the hits (each hit corresponds to a fiber)
  for (const auto& scintHit : simHits) {

    if (!scintHit.isAvailable()) {
      error() << "Hit not available!" << endmsg;
      continue;
    }

    std::vector<double> vecTimes;
    std::vector<double> vecWavelens;

    // loop through the hit contributions
    for (auto contrib = scintHit.contributions_begin(); contrib != scintHit.contributions_end(); ++contrib) {

      if (!contrib->isAvailable()) {
        error() << "Hit contribution not available!" << endmsg;
        continue;
      }
      const double npe = contrib->getEnergy(); // in IDEA_o2 this is the number of photo-electrons

      vecTimes.reserve(vecTimes.size() + npe); // resize vector to store new p.e.
      vecWavelens.reserve(vecTimes.size() + npe);

      // get the photon arrival time at SiPM
      const double time = contrib->getTime(); // in IDEA_o2 this it the toa at SiPM in ns
      double thisTime = time;
      // if it is a scintillation hit, add scintillation decay time
      if (!m_isCherenkov) {
        thisTime += m_engine.Exp(m_scintDecaytime.value());
      }

      for (unsigned int pe = 0; pe < npe; pe++) {
        vecTimes.push_back(thisTime);
        vecWavelens.push_back(1.); // photons wavelengths are dummy, we are using photo-electrons
      }
    } // contrib

    // sort time of arrivals from shortest to longest
    std::sort(vecTimes.begin(), vecTimes.end());

    m_sensor->resetState();
    m_sensor->addPhotons(vecTimes, vecWavelens); // Sets photon times & wavelengths
    m_sensor->runEvent();                        // Runs the simulation

    // Using only analog signal (ADC conversion is still experimental)
    const sipm::SiPMAnalogSignal anaSignal = m_sensor->signal();

    // if the signal never exceeds the threshold, it will return -1.
    const double integral =
        std::max(0., anaSignal.integral(m_gateStart, m_gateL, m_thres));           // (intStart, intGate, threshold)
    const double toa = std::max(0., anaSignal.toa(m_gateStart, m_gateL, m_thres)); // (intStart, intGate, threshold)

    // Apply threshold to save digi hit.
    // By default this is a zero suppression cut, can the changed to apply 1-suppression, 2-suppression, ...
    // Note: check the sipm gain applied in SimSipm, by default it is 1.
    if (integral < m_suppressionThreshold)
      continue;

    // NOTE: I noticed very few digi hits have TOA > ms, clearly a bug to be investigated
    if (toa > 1e6)
      continue;

    auto digiHit = digiHits.create();
    digiHit.setEnergy(integral /* * m_scaleADC.value()*/); // convertition from ADC to GeV can be added here
    digiHit.setEnergyError(/*m_scaleADC.value() **/ std::sqrt(integral));
    digiHit.setPosition(scintHit.getPosition());
    digiHit.setCellID(scintHit.getCellID());
    // Toa and m_gateStart are in ns
    digiHit.setTime(toa + m_gateStart);

    auto link = links.create();
    link.setFrom(digiHit);
    link.setTo(scintHit);
  }

  return std::make_tuple(std::move(digiHits), std::move(links));
}
