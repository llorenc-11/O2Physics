// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file jetHadronsPid.cxx
/// \brief Analysis of hadrons in jets
/// \author Leonard Lorenc, WUT Warsaw, leonard.lorenc@cern.ch
/// \author Aleksandra Mulewicz, WUT Warsaw, aleksandra.mulewicz@cern.ch
/// \author Hubert Zalewski, WUT Warsaw, hubert.zalewski@cern.ch
/// \author Janik Małgorzata malgorzata.janik@cern.ch
/// \author Daniela Ruggiano daniela.ruggiano@cern.ch

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>
#include <TVector3.h>

#include <Rtypes.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;

using StandardEvents = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;

using JetEvents = soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionPIs>;

using MCDJetEvents = soa::Join<aod::JetCollisions, aod::JCollisionPIs>;

using ChargedMCDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;

using HadronTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA,
                               aod::pidTPCFullPi, aod::pidTOFFullPi,
                               aod::pidTPCFullKa, aod::pidTOFFullKa,
                               aod::pidTPCFullPr, aod::pidTOFFullPr,
                               aod::pidTOFbeta>;

using HadronTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA,
                                 aod::pidTPCFullPi, aod::pidTOFFullPi,
                                 aod::pidTPCFullKa, aod::pidTOFFullKa,
                                 aod::pidTPCFullPr, aod::pidTOFFullPr,
                                 aod::pidTOFbeta, aod::McTrackLabels>;

struct JetHadronsPid {

  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<std::string> eventSelections{"eventSelections", "sel8,vtxITSTPC,noITSROF,noTF,noPileup,goodZvtx", "choose event selection"};
  std::vector<int> eventSelectionBits;

  Configurable<bool> isppRefAnalysis{"isppRefAnalysis", true, "Is ppRef analysis"};
  Configurable<double> cfgEtaJetMax{"cfgEtaJetMax", 0.5, "max jet eta"};

  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "Reject events near the ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "Reject events near the TF border"};
  Configurable<bool> requireVtxITSTPC{"requireVtxITSTPC", true, "Require at least one ITS-TPC matched track"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "Reject events with same-bunch pileup collisions"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "Require consistent FT0 vs PV z-vertex"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "Require at least one vertex track matched to TOF"};

  Configurable<double> minJetPt{"minJetPt", 10, "Minimum pt of the jet after bkg subtraction"};
  Configurable<double> maxJetPt{"maxJetPt", 1e+06, "Maximum pt of the jet after bkg subtraction"};
  Configurable<double> rJet{"rJet", 0.4, "Jet resolution parameter R"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<bool> applyAreaCut{"applyAreaCut", false, "apply area cut"};
  Configurable<double> minNormalizedJetArea{"minNormalizedJetArea", 0.6, "Minimum normalized area cut to reject fake jets"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from the edge"};

  Configurable<bool> requirePvContributor{"requirePvContributor", false, "require that the track is a PV contributor"};
  Configurable<int> minItsNclusters{"minItsNclusters", 5, "minimum number of ITS clusters"};
  Configurable<int> minTpcNcrossedRows{"minTpcNcrossedRows", 70, "minimum number of TPC crossed pad rows"};
  Configurable<double> minChiSquareTpc{"minChiSquareTpc", 0.0, "minimum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareTpc{"maxChiSquareTpc", 4.0, "maximum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareIts{"maxChiSquareIts", 36.0, "maximum ITS chi^2/Ncls"};
  Configurable<double> minPt{"minPt", 0.05, "minimum pt of the tracks"};
  Configurable<double> maxPt{"maxPt", 4.0, "maximum pt of the tracks for PID analysis"};
  Configurable<double> minEta{"minEta", -0.8, "minimum eta"};
  Configurable<double> maxEta{"maxEta", +0.8, "maximum eta"};
  Configurable<double> maxDcaxy{"maxDcaxy", 0.2, "Maximum DCAxy"};
  Configurable<double> maxDcaz{"maxDcaz", 0.1, "Maximum DCAz"};
  Configurable<bool> setMCDefaultItsParams{"setMCDefaultItsParams", true, "set MC default parameters"};

  o2::aod::ITSResponse itsResponse;

  Preslice<HadronTracks> tracksPerCollision = aod::track::collisionId;
  Preslice<HadronTracksMC> mcTracksPerCollision = aod::track::collisionId;
  Preslice<aod::McParticles> mcParticlesPerCollision = aod::mcparticle::mcCollisionId;

  struct : ConfigurableGroup {
    Configurable<int> pidMethod{"pidMethod", 2, "Variable for choosing pid method"};
    Configurable<float> rejectionSigma{"rejectionSigma", 3.0, "Rejection sigma for tof and tpc pid"};
    Configurable<float> ptThresholdPion{"ptThresholdPion", 0.75, "Threshold for pt for different pid methods"};
    Configurable<float> ptThresholdKaon{"ptThresholdKaon", 0.5, "Threshold for pt for different pid methods"};
    Configurable<float> ptThresholdProton{"ptThresholdProton", 0.75, "Threshold for pt for different pid methods"};
    Configurable<float> nSigmaCut{"nSigmaCut", 3.0, "Acceptence cut for pid methods"};
    Configurable<float> minPtPion{"minPtPion", 0.2, "Minimal Pion pt"};
    Configurable<float> maxPtPion{"maxPtPion", 4.0, "Maximum Pion pt"};
    Configurable<float> minPtKaon{"minPtKaon", 0.3, "Minimal Kaon pt"};
    Configurable<float> maxPtKaon{"maxPtKaon", 3.0, "Maximum Kaon pt"};
    Configurable<float> minPtProton{"minPtProton", 0.5, "Minimal Proton pt"};
    Configurable<float> maxPtProton{"maxPtProton", 4.0, "Maximum Proton pt"};
  } cfg;

  void init(InitContext const&)
  {
    if (setMCDefaultItsParams) {
      itsResponse.setMCDefaultParameters();
    }

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    
    //////////////////////////////////////////////
    //                AXIS SPEC
    //////////////////////////////////////////////

    AxisSpec axisMomentumTPCSignal{300, 0.05, 30.0, "#it{p} (GeV/#it{c})"};
    AxisSpec axisTPCSignal{400, 20.0, 2000.0, "TPC d#it{E}/d#it{x} signal (a.u.)"};
    AxisSpec axisMomentumTOF{300, 0.05, 30.0, "#it{p} (GeV/#it{c})"};
    AxisSpec axisTOFBeta{240, 0.0, 1.2, "#beta_{TOF}"};
    AxisSpec axisDeltaEtaCone{80, -rJet, rJet, "#Delta#eta = #eta_{track} - #eta_{jet}"};
    AxisSpec axisDeltaPhiCone{80, -rJet, rJet, "#Delta#varphi = #varphi_{track} - #varphi_{jet}"};
    AxisSpec axisPtTrackCone{100, 0.0, 10.0, "#it{p}_{T}^{track} (GeV/#it{c})"};
    AxisSpec axisMultiplicity{31, -0.5, 30.5, "N per jet"};
    AxisSpec axisNPions{21, -0.5, 20.5, "N_{#pi} per jet"};
    AxisSpec axisNKaons{11, -0.5, 10.5, "N_{K} per jet"};
    AxisSpec axisNProtons{11, -0.5, 10.5, "N_{p} per jet"};
    AxisSpec axisMomentumMC{300, 0.05, 30.0, "#it{p} (GeV/#it{c})"};
    AxisSpec axisTPCSignalMC{400, 20.0, 2000.0, "TPC d#it{E}/d#it{x} signal (a.u.)"};
    AxisSpec axisTOFBetaMC{240, 0.0, 1.2, "#beta_{TOF}"};
    AxisSpec axisMCDDeltaEtaCone{80, -static_cast<double>(rJet), static_cast<double>(rJet), "#Delta#eta = #eta_{track} - #eta_{axis}"};
    AxisSpec axisMCDDeltaPhiCone{80, -static_cast<double>(rJet), static_cast<double>(rJet), "#Delta#varphi = #varphi_{track} - #varphi_{axis}"};
    AxisSpec axisMCDPtTrackCone{100, 0.0, 10.0, "#it{p}_{T}^{track} (GeV/#it{c})"};
    AxisSpec axisMCDMultiplicity{31, -0.5, 30.5, "N per jet"};
    AxisSpec axisMCDNPions{21, -0.5, 20.5, "N_{#pi} per jet"};
    AxisSpec axisMCDNKaons{11, -0.5, 10.5, "N_{K} per jet"};
    AxisSpec axisMCDNProtons{11, -0.5, 10.5, "N_{p} per jet"};
    AxisSpec axisTruthPt{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisTruthEta{80, -0.8, 0.8, "#eta"};
    AxisSpec axisTruthPhi{72, 0.0, TwoPI, "#varphi"};
    AxisSpec axisTruthMultiplicity{201, -0.5, 200.5, "N_{particles}"};

    //////////////////////////////////////////////
    //                   PURE
    //////////////////////////////////////////////

    registryData.add("data/pure/tpc_signal_vs_p", "TPC n#sigma_{#pi} vs transverse momentum ;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axisMomentumTPCSignal, axisTPCSignal});
    registryData.add("data/pure/pions/pion_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected pion candidates;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {{300, 0.05, 30.0, "#it{p} (GeV/#it{c})"}, {400, 20.0, 2000.0, "TPC d#it{E}/d#it{x} signal (a.u.)"}});
    registryData.add("data/pure/kaons/kaon_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected kaon candidates;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F,{{300, 0.05, 30.0, "#it{p} (GeV/#it{c})"}, {400, 20.0, 2000.0, "TPC d#it{E}/d#it{x} signal (a.u.)"}});
    registryData.add("data/pure/protons/proton_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected proton candidates;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {{300, 0.05, 30.0, "#it{p} (GeV/#it{c})"}, {400, 20.0, 2000.0, "TPC d#it{E}/d#it{x} signal (a.u.)"}});

    registryData.add("data/pure/tof_beta_vs_p", "TOF #beta vs momentum;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axisMomentumTOF, axisTOFBeta});
    registryData.add("data/pure/pions/pion_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected pion candidates;#it{p} (GeV/#it{c});#beta_{TOF}" "#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axisMomentumTOF, axisTOFBeta});
    registryData.add("data/pure/kaons/kaon_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected kaon candidates;#it{p} (GeV/#it{c});#beta_{TOF}",  HistType::kTH2F, {axisMomentumTOF, axisTOFBeta});
    registryData.add("data/pure/protons/proton_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected proton candidates;#it{p} (GeV/#it{c});#beta_{TOF}" "#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axisMomentumTOF, axisTOFBeta});

    // Pions
    registryData.add("data/pure/pions/pion_pure_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected pion candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{#pi}"}});
    registryData.add("data/pure/pions/pion_pure_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected pion candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}^{#pi}"}});
    registryData.add("data/pure/pions/pion_pure_pt", "Pion pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/pure/pions/pion_pure_eta", "Pion Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/pure/pions/pion_pure_dcaxy", "Pion DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/pure/pions/pion_pure_dcaz", "Pion DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/pure/pions/pion_pure_eta_phi", "Azimuthal angle vs pseudorapidity for PID-selected pion candidates", HistType::kTH2F, {{72, 0.0, TwoPI}, {80, -0.8, 0.8}});

    registryData.add("data/pure/pions/pos/pion_pure_pos_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{+} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{#pi}"}});
    registryData.add("data/pure/pions/pos/pion_pure_pos_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{+} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}^{#pi}"}});
    registryData.add("data/pure/pions/pos/pion_pure_pos_pt", "#pi^{+} pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/pure/pions/pos/pion_pure_pos_eta", "#pi^{+} Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/pure/pions/pos/pion_pure_pos_dcaxy", "#pi^{+} DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/pure/pions/pos/pion_pure_pos_dcaz", "#pi^{+} DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/pure/pions/pos/pion_pure_pos_phi_vs_pt", "Pseudorapidity vs azimuthal angle for PID-selected #pi^{+} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});
    registryData.add("data/pure/pions/pos/pion_pure_pos_eta_phi", "Azimuthal angle vs transverse momentum for PID-selected #pi^{+} candidates", HistType::kTH2F, {{72, 0.0, TwoPI, "#varphi"}, {80, -0.8, 0.8, "#eta"}});

    registryData.add("data/pure/pions/neg/pion_pure_neg_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{-} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{#pi}"}});
    registryData.add("data/pure/pions/neg/pion_pure_neg_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{-} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}^{#pi}"}});
    registryData.add("data/pure/pions/neg/pion_pure_neg_pt", "#pi^{-} pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/pure/pions/neg/pion_pure_neg_eta", "#pi^{-} Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/pure/pions/neg/pion_pure_neg_dcaxy", "#pi^{-} DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/pure/pions/neg/pion_pure_neg_dcaz", "#pi^{-} DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/pure/pions/neg/pion_pure_neg_phi_vs_pt", "Pseudorapidity vs azimuthal angle for PID-selected #pi^{-} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});
    registryData.add("data/pure/pions/neg/pion_pure_neg_eta_phi", "Azimuthal angle vs transverse momentum for PID-selected #pi^{-} candidates", HistType::kTH2F, {{72, 0.0, TwoPI, "#varphi"}, {80, -0.8, 0.8, "#eta"}});

    // Kaons
    registryData.add("data/pure/kaons/kaon_pure_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected kaon candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{K}"}});
    registryData.add("data/pure/kaons/kaon_pure_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected kaon candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/pure/kaons/kaon_pure_pt", "Kaon pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/pure/kaons/kaon_pure_eta", "Kaon Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/pure/kaons/kaon_pure_dcaz", "Kaon DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/pure/kaons/kaon_pure_dcaxy", "Kaon DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/pure/kaons/kaon_pure_eta_phi", "Azimuthal angle vs transverse momentum for PID-selected Kaon candidates", HistType::kTH2F, {{72, 0.0, TwoPI}, {80, -0.8, 0.8}});

    registryData.add("data/pure/kaons/pos/kaon_pure_pos_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected K^{+} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{K}"}});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected K^{+} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_pt", "K^{+} pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_eta", "K^{+} Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_dcaxy", "K^{+} DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_dcaz", "K^{+} DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_phi_vs_pt", "Pseudorapidity vs azimuthal angle for PID-selected K^{+} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_eta_phi", "Azimuthal angle vs transverse momentum for PID-selected K^{+} candidates", HistType::kTH2F, {{72, 0.0, TwoPI, "#varphi"}, {80, -0.8, 0.8, "#eta"}});

    registryData.add("data/pure/kaons/neg/kaon_pure_neg_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected K^{-} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{K}"}});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected K^{-} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_pt", "K^{-} pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_eta", "K^{-} Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_dcaxy", "K^{-} DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_dcaz", "K^{-} DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_phi_vs_pt", "Pseudorapidity vs azimuthal angle for PID-selected K^{-} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_eta_phi", "Azimuthal angle vs transverse momentum for PID-selected K^{-} candidates", HistType::kTH2F, {{72, 0.0, TwoPI, "#varphi"}, {80, -0.8, 0.8, "#eta"}});

    // Protons
    registryData.add("data/pure/protons/proton_pure_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected protons + antiprotons candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/pure/protons/proton_pure_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected protons + antiprotons candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/pure/protons/proton_pure_pt", "Proton pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/pure/protons/proton_pure_eta", "Proton Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/pure/protons/proton_pure_dcaz", "Proton DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/pure/protons/proton_pure_dcaxy", "Proton DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/pure/protons/proton_pure_eta_phi", "Azimuthal angle vs transverse momentum for PID-selected protons + antiprotons candidates", HistType::kTH2F, {{72, 0.0, TwoPI}, {80, -0.8, 0.8}});

    registryData.add("data/pure/protons/pos/proton_pure_pos_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected p candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/pure/protons/pos/proton_pure_pos_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected p candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/pure/protons/pos/proton_pure_pos_pt", "p pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/pure/protons/pos/proton_pure_pos_eta", "p Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/pure/protons/pos/proton_pure_pos_dcaxy", "p DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/pure/protons/pos/proton_pure_pos_dcaz", "p DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/pure/protons/pos/proton_pure_pos_phi_vs_pt", "Pseudorapidity vs azimuthal angle for PID-selected p candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});
    registryData.add("data/pure/protons/pos/proton_pure_pos_eta_phi", "Azimuthal angle vs transverse momentum for PID-selected p candidates", HistType::kTH2F, {{72, 0.0, TwoPI, "#varphi"}, {80, -0.8, 0.8, "#eta"}});

    registryData.add("data/pure/protons/neg/proton_pure_neg_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected #bar{p} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/pure/protons/neg/proton_pure_neg_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected #bar{p} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/pure/protons/neg/proton_pure_neg_pt", "#bar{p} pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/pure/protons/neg/proton_pure_neg_eta", "#bar{p} Eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/pure/protons/neg/proton_pure_neg_dcaxy", "#bar{p} DCAxy", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/pure/protons/neg/proton_pure_neg_dcaz", "#bar{p} DCAz", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/pure/protons/neg/proton_pure_neg_phi_vs_pt", "Pseudorapidity vs azimuthal angle for PID-selected #bar{p} candidates", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});
    registryData.add("data/pure/protons/neg/proton_pure_neg_eta_phi", "Azimuthal angle vs transverse momentum for PID-selected #bar{p} candidates", HistType::kTH2F, {{72, 0.0, TwoPI, "#varphi"}, {80, -0.8, 0.8, "#eta"}});

    //////////////////////////////////////////////
    //                COLLISIONS
    //////////////////////////////////////////////

    registryData.add("data/pure/collisions/z_vertex", "Collision vertex;z_{vtx} (cm);N_{collisions}", HistType::kTH1F, {{200, -20.0, 20.0, "z_{vtx} (cm)"}});
    registryData.add("data/pure/collisions/x_vertex", "Collision vertex;x_{vtx} (cm);N_{collisions}", HistType::kTH1F, {{200, -1.0, 1.0, "x_{vtx} (cm)"}});
    registryData.add("data/pure/collisions/y_vertex", "Collision vertex;y_{vtx} (cm);N_{collisions}", HistType::kTH1F, {{200, -1.0, 1.0, "y_{vtx} (cm)"}});
    registryData.add("data/pure/collisions/xy_vertex", "Collision vertex;x_{vtx} (cm);y_{vtx} (cm)", HistType::kTH2F, {{200, -1.0, 1.0, "x_{vtx} (cm)"}, {200, -1.0, 1.0, "y_{vtx} (cm)"}});
    registryData.add("data/pure/collisions/n_contributors", "PV contributors;N_{contributors};N_{collisions}", HistType::kTH1I, {{201, -0.5, 200.5, "N_{contributors}"}});

    //////////////////////////////////////////////
    //                   JETS
    //////////////////////////////////////////////

    registryData.add("data/jets/collision_multiplicity", "Number of Jets in Collison;Number of Jets;Number of Collisions", HistType::kTH1F, {{30, 0.0, 30.5}});
    registryData.add("data/jets/n_events", "Event counter", HistType::kTH1F, {{1, 0.5, 1.5, "N_{events}"}});
    registryData.add("data/jets/n_events_raw", "All events", HistType::kTH1F, {{1, 0.5, 1.5, ""}});

    registryData.add("data/jets/jet_pt", "Jet pT ", HistType::kTH1F, {{100, 0.0, 20.0, "#it{p}_{T}^{raw} (GeV/#it{c})"}});
    registryData.add("data/jets/jet_pt_subtracted", "Jet pT subtracted", HistType::kTH1F, {{200, 0.0, 10.0, "#it{p}_{T}^{sub} (GeV/#it{c})"}});
    registryData.add("data/jets/jet_pt_raw_vs_sub", "Raw vs sub jet pT", HistType::kTH2F, {{200, 0, 200}, {200, 0, 200}});
    registryData.add("data/jets/jet_eta", "Jet eta", HistType::kTH1F, {{100, -1.0, 1.0, "#eta_{jet}"}});
    registryData.add("data/jets/jet_phi", "Jet phi", HistType::kTH1F, {{100, 0.0, TwoPI, "#phi_{jet}"}});
    registryData.add("data/jets/jet_area", "Jet area", HistType::kTH1F, {{100, 0.0, 1.5, "Area"}});
    registryData.add("data/jets/jet_n_constituents", "Jet multiplicity", HistType::kTH1I, {{100, 0, 30, "N_{constituents}"}});
    registryData.add("data/jets/z_vtx", "Z-Vertex Distribution", HistType::kTH1F, {{200, -20.0, 20.0, "Z-Vertex (cm)"}});
    registryData.add("data/jets/tpc_signal_vs_p", "TPC signal for tracks in selected jets;" "#it{p} (GeV/#it{c});" "TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axisMomentumTPCSignal, axisTPCSignal});
    registryData.add("data/jets/tof_beta_vs_p", "TOF #beta for tracks in selected jets;" "#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axisMomentumTOF, axisTOFBeta});
    
   // Jet Cone

    registryData.add("data/jets/cone/all_tracks_2d", "Track distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axisDeltaEtaCone, axisDeltaPhiCone});
    registryData.add("data/jets/cone/pions/pion_pos_2d", "PID-selected #pi^{+} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axisDeltaEtaCone, axisDeltaPhiCone});
    registryData.add("data/jets/cone/pions/pion_neg_2d", "PID-selected #pi^{-} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axisDeltaEtaCone, axisDeltaPhiCone});
    registryData.add("data/jets/cone/kaons/kaon_pos_2d", "PID-selected K^{+} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axisDeltaEtaCone, axisDeltaPhiCone});
    registryData.add("data/jets/cone/kaons/kaon_neg_2d", "PID-selected K^{-} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axisDeltaEtaCone, axisDeltaPhiCone});
    registryData.add("data/jets/cone/protons/proton_2d", "PID-selected proton distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axisDeltaEtaCone, axisDeltaPhiCone});
    registryData.add("data/jets/cone/protons/antiproton_2d", "PID-selected antiproton distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axisDeltaEtaCone, axisDeltaPhiCone});

    registryData.add("data/jets/cone/all_tracks", "Track distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axisDeltaEtaCone, axisDeltaPhiCone, axisPtTrackCone});
    registryData.add("data/jets/cone/pions/pion_pos", "PID-selected #pi^{+} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axisDeltaEtaCone, axisDeltaPhiCone, axisPtTrackCone});
    registryData.add("data/jets/cone/pions/pion_neg", "PID-selected #pi^{-} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axisDeltaEtaCone, axisDeltaPhiCone, axisPtTrackCone});
    registryData.add("data/jets/cone/kaons/kaon_pos", "PID-selected K^{+} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axisDeltaEtaCone, axisDeltaPhiCone, axisPtTrackCone});
    registryData.add("data/jets/cone/kaons/kaon_neg", "PID-selected K^{-} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axisDeltaEtaCone, axisDeltaPhiCone, axisPtTrackCone});
    registryData.add("data/jets/cone/protons/proton", "PID-selected proton distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axisDeltaEtaCone, axisDeltaPhiCone, axisPtTrackCone});
    registryData.add("data/jets/cone/protons/antiproton", "PID-selected antiproton distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axisDeltaEtaCone, axisDeltaPhiCone, axisPtTrackCone});

    registryData.add("data/jets/pions/pos/n_pions_per_jet", "Multiplicity of PID-selected #pi^{+} candidates in selected jets;N_{#pi^{+}};N_{jets}", HistType::kTH1I, {axisMultiplicity});
    registryData.add("data/jets/kaons/pos/n_kaons_per_jet", "Multiplicity of PID-selected K^{+} candidates in selected jets;N_{K^{+}};N_{jets}", HistType::kTH1I, {axisMultiplicity});
    registryData.add("data/jets/protons/pos/n_protons_per_jet", "Multiplicity of PID-selected protons in selected jets;N_{p};N_{jets}", HistType::kTH1I, {axisMultiplicity});
    registryData.add("data/jets/pions/neg/n_pions_per_jet", "Multiplicity of PID-selected #pi^{-} candidates in selected jets;N_{#pi^{-}};N_{jets}", HistType::kTH1I, {axisMultiplicity});
    registryData.add("data/jets/kaons/neg/n_kaons_per_jet", "Multiplicity of PID-selected K^{-} candidates in selected jets;N_{K^{-}};N_{jets}", HistType::kTH1I, {axisMultiplicity});
    registryData.add("data/jets/protons/neg/n_protons_per_jet", "Multiplicity of PID-selected antiprotons in selected jets;N_{#bar{p}};N_{jets}", HistType::kTH1I, {axisMultiplicity});
    registryData.add("data/jets/composition/particle_multiplicity_per_jet", "Identified-particle multiplicity in selected jets;Particle species;N_{particles} per jet", HistType::kTH2I, {{3, 0.5, 3.5, "Particle species"}, {31, -0.5, 30.5, "N_{particles} per jet"}});

    auto hComposition = registryData.get<TH2>(HIST("data/jets/composition/particle_multiplicity_per_jet"));
    hComposition->GetXaxis()->SetBinLabel(1, "#pi");
    hComposition->GetXaxis()->SetBinLabel(2, "K");
    hComposition->GetXaxis()->SetBinLabel(3, "p");
    hComposition->SetOption("COLZ");
    hComposition->SetStats(false);

    registryData.add("data/jets/composition/nKaons_vs_nProtons", "Particle multiplicity correlation;N_{K};N_{p};N_{jets}", HistType::kTH2F, {axisNKaons, axisNProtons});
    registryData.add("data/jets/composition/nPions_vs_nKaons", "Particle multiplicity correlation;N_{#pi};N_{K};N_{jets}", HistType::kTH2F, {axisNPions, axisNKaons});
    registryData.add("data/jets/composition/nPions_vs_nProtons", "Particle multiplicity correlation;N_{#pi};N_{p};N_{jets}", HistType::kTH2F, {axisNPions, axisNProtons});
    registryData.add("data/jets/composition/nPions_nKaons_nProtons", "Jet particle composition;N_{#pi};N_{K};N_{p}", HistType::kTH3F, {axisNPions, axisNKaons, axisNProtons});

    registryData.add("data/jets/composition/nKaons_vs_nPions_1proton", "Kaons vs pions in jets with exactly 1 p/#bar{p};N_{K};N_{#pi};N_{jets}", HistType::kTH2I, {axisNKaons, axisNPions});
    registryData.add("data/jets/composition/nKaons_vs_nPions_2protons", "Kaons vs pions in jets with exactly 2 p/#bar{p};N_{K};N_{#pi};N_{jets}", HistType::kTH2I, {axisNKaons, axisNPions});
    registryData.add("data/jets/composition/nKaons_vs_nPions_3protons", "Kaons vs pions in jets with exactly 3 p/#bar{p};N_{K};N_{#pi};N_{jets}", HistType::kTH2I, {axisNKaons, axisNPions});
    registryData.add("data/jets/composition/nKaons_vs_nPions_4protons", "Kaons vs pions in jets with exactly 4 p/#bar{p};N_{K};N_{#pi};N_{jets}", HistType::kTH2I, {axisNKaons, axisNPions});

    // Pions

    registryData.add("data/jets/pions/pion_tpc_signal_vs_p", "TPC signal after pion PID, tracks in jets;" "#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {{300, 0.05, 30.0, "#it{p} (GeV/#it{c})"}, {400, 20.0, 2000.0, "TPC d#it{E}/d#it{x} signal (a.u.)"}});
    registryData.add("data/jets/pions/pion_tof_beta_vs_p", "TOF #beta after pion PID, tracks in jets;" "#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axisMomentumTOF, axisTOFBeta});

    registryData.add("data/jets/pions/pos/pion_jet_pos_tpc", "TPC #pi^{+} PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/jets/pions/pos/pion_jet_pos_tof", "TOF #pi^{+} PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/jets/pions/pos/pion_jet_pos_pt", "#pi^{+} pT in Jets", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/jets/pions/pos/pion_jet_pos_eta", "#pi^{+} Eta in Jets", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/jets/pions/pos/pion_jet_pos_dcaxy", "#pi^{+} DCAxy in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/jets/pions/pos/pion_jet_pos_dcaz", "#pi^{+} DCAz in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/jets/pions/pos/pion_jet_pos_phi_vs_pt", "#pi^{+} in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});

    registryData.add("data/jets/pions/neg/pion_jet_neg_tpc", "TPC #pi^{-} PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/jets/pions/neg/pion_jet_neg_tof", "TOF #pi^{-} PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/jets/pions/neg/pion_jet_neg_pt", "#pi^{-} pT in Jets", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/jets/pions/neg/pion_jet_neg_eta", "#pi^{-} Eta in Jets", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/jets/pions/neg/pion_jet_neg_dcaxy", "#pi^{-} DCAxy in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/jets/pions/neg/pion_jet_neg_dcaz", "#pi^{-} DCAz in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/jets/pions/neg/pion_jet_neg_phi_vs_pt", "#pi^{-} in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});

    // Kaons

    registryData.add("data/jets/kaons/kaon_tpc_signal_vs_p", "TPC signal after kaon PID, tracks in jets;" "#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {{300, 0.05, 30.0, "#it{p} (GeV/#it{c})"}, {400, 20.0, 2000.0, "TPC d#it{E}/d#it{x} signal (a.u.)"}});
    registryData.add("data/jets/kaons/kaon_tof_beta_vs_p", "TOF #beta after kaon PID, tracks in jets;" "#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axisMomentumTOF, axisTOFBeta});

    registryData.add("data/jets/kaons/pos/kaon_jet_pos_tpc", "TPC K^{+} PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_tof", "TOF K^{+} PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_pt", "K^{+} pT in Jets", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_eta", "K^{+} Eta in Jets", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_dcaxy", "K^{+} DCAxy in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_dcaz", "K^{+} DCAz in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_phi_vs_pt", "K^{+} in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});
    
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_tpc", "TPC K^{-} PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_tof", "TOF K^{-} PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_pt", "K^{-} pT in Jets", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_eta", "K^{-} Eta in Jets", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_dcaxy", "K^{-} DCAxy in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_dcaz", "K^{-} DCAz in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_phi_vs_pt", "K^{-} in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});

    // Protons
    registryData.add("data/jets/protons/proton_tpc_signal_vs_p", "TPC signal after proton PID, tracks in jets;" "#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {{300, 0.05, 30.0, "#it{p} (GeV/#it{c})"}, {400, 20.0, 2000.0, "TPC d#it{E}/d#it{x} signal (a.u.)"}});
    registryData.add("data/jets/protons/proton_tof_beta_vs_p", "TOF #beta after proton PID, tracks in jets;" "#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axisMomentumTOF, axisTOFBeta});

    registryData.add("data/jets/protons/pos/proton_jet_pos_tpc", "TPC p PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/jets/protons/pos/proton_jet_pos_tof", "TOF p PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/jets/protons/pos/proton_jet_pos_pt", "p pT in Jets", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/jets/protons/pos/proton_jet_pos_eta", "p Eta in Jets", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/jets/protons/pos/proton_jet_pos_dcaxy", "p DCAxy in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/jets/protons/pos/proton_jet_pos_dcaz", "p DCAz in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/jets/protons/pos/proton_jet_pos_phi_vs_pt", "p in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});

    registryData.add("data/jets/protons/neg/proton_jet_neg_tpc", "TPC #bar{p} PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/jets/protons/neg/proton_jet_neg_tof", "TOF #bar{p} PID in Jets", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/jets/protons/neg/proton_jet_neg_pt", "#bar{p} pT in Jets", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/jets/protons/neg/proton_jet_neg_eta", "#bar{p} Eta in Jets", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/jets/protons/neg/proton_jet_neg_dcaxy", "#bar{p} DCAxy in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/jets/protons/neg/proton_jet_neg_dcaz", "#bar{p} DCAz in Jets", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/jets/protons/neg/proton_jet_neg_phi_vs_pt", "p in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});

    //////////////////////////////////////////////
    //                    UE
    //////////////////////////////////////////////

    // Pions

    registryData.add("data/ue/pions/pion_ue_tpc", "TPC Pion PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/ue/pions/pion_ue_tof", "TOF Pion PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/ue/pions/pion_ue_pt", "Pion pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/ue/pions/pion_ue_eta", "Pion Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/ue/pions/pion_ue_dcaxy", "Pion DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/ue/pions/pion_ue_dcaz", "Pion DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("data/ue/pions/pos/pion_ue_pos_tpc", "TPC #pi^{+} PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/ue/pions/pos/pion_ue_pos_tof", "TOF #pi^{+} PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/ue/pions/pos/pion_ue_pos_pt", "#pi^{+} pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/ue/pions/pos/pion_ue_pos_eta", "#pi^{+} Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/ue/pions/pos/pion_ue_pos_dcaxy", "#pi^{+} DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/ue/pions/pos/pion_ue_pos_dcaz", "#pi^{+} DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("data/ue/pions/neg/pion_ue_neg_tpc", "TPC #pi^{-} PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/ue/pions/neg/pion_ue_neg_tof", "TOF #pi^{-} PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/ue/pions/neg/pion_ue_neg_pt", "#pi^{-} pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/ue/pions/neg/pion_ue_neg_eta", "#pi^{-} Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/ue/pions/neg/pion_ue_neg_dcaxy", "#pi^{-} DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/ue/pions/neg/pion_ue_neg_dcaz", "#pi^{-} DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    // Kaons

    registryData.add("data/ue/kaons/kaon_ue_tpc", "TPC Kaon PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/ue/kaons/kaon_ue_tof", "TOF Kaon PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/ue/kaons/kaon_ue_pt", "Kaon pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/ue/kaons/kaon_ue_eta", "Kaon Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/ue/kaons/kaon_ue_dcaxy", "Kaon DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/ue/kaons/kaon_ue_dcaz", "Kaon DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("data/ue/kaons/pos/kaon_ue_pos_tpc", "TPC K^{+} PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/ue/kaons/pos/kaon_ue_pos_tof", "TOF K^{+} PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/ue/kaons/pos/kaon_ue_pos_pt", "K^{+} pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/ue/kaons/pos/kaon_ue_pos_eta", "K^{+} Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/ue/kaons/pos/kaon_ue_pos_dcaxy", "K^{+} DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/ue/kaons/pos/kaon_ue_pos_dcaz", "K^{+} DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("data/ue/kaons/neg/kaon_ue_neg_tpc", "TPC K^{-} PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/ue/kaons/neg/kaon_ue_neg_tof", "TOF K^{-} PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/ue/kaons/neg/kaon_ue_neg_pt", "K^{-} pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/ue/kaons/neg/kaon_ue_neg_eta", "K^{-} Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/ue/kaons/neg/kaon_ue_neg_dcaxy", "K^{-} DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/ue/kaons/neg/kaon_ue_neg_dcaz", "K^{-} DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    // Protons

    registryData.add("data/ue/protons/proton_ue_tpc", "TPC Proton PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/ue/protons/proton_ue_tof", "TOF Proton PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/ue/protons/proton_ue_pt", "Proton pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/ue/protons/proton_ue_eta", "Proton Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/ue/protons/proton_ue_dcaxy", "Proton DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/ue/protons/proton_ue_dcaz", "Proton DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("data/ue/protons/pos/proton_ue_pos_tpc", "TPC p PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/ue/protons/pos/proton_ue_pos_tof", "TOF p PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/ue/protons/pos/proton_ue_pos_pt", "p pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/ue/protons/pos/proton_ue_pos_eta", "p Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/ue/protons/pos/proton_ue_pos_dcaxy", "p DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/ue/protons/pos/proton_ue_pos_dcaz", "p DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});

    registryData.add("data/ue/protons/neg/proton_ue_neg_tpc", "TPC #bar{p} PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}"}});
    registryData.add("data/ue/protons/neg/proton_ue_neg_tof", "TOF #bar{p} PID in UE", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}"}});
    registryData.add("data/ue/protons/neg/proton_ue_neg_pt", "#bar{p} pT in UE", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("data/ue/protons/neg/proton_ue_neg_eta", "#bar{p} Eta in UE", HistType::kTH1F, {{100, -1.0, 1.0, "#eta"}});
    registryData.add("data/ue/protons/neg/proton_ue_neg_dcaxy", "#bar{p} DCAxy in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    registryData.add("data/ue/protons/neg/proton_ue_neg_dcaz", "#bar{p} DCAz in UE", HistType::kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    registryData.add("data/ue/tracks_n_in_ue", "Number of tracks in UE", HistType::kTH1I, {{100, 0, 100, "N_{tracks}"}});

    //////////////////////////////////////////////
    //             MC RECONSTRUCTION
    //////////////////////////////////////////////

    registryData.add("mc/reconstruction/collisions/n_events", "Selected reconstructed MC collisions;event;N_{collisions}", HistType::kTH1I, {{1, 0.5, 1.5, "event"}});
    registryData.add("mc/reconstruction/collisions/z_vertex", "Reconstructed MC collision vertex;z_{vtx}^{rec} (cm);N_{collisions}", HistType::kTH1F, {{200, -20.0, 20.0, "z_{vtx}^{rec} (cm)"}});
    registryData.add("mc/reconstruction/collisions/x_vertex", "Reconstructed MC collision vertex;x_{vtx}^{rec} (cm);N_{collisions}", HistType::kTH1F, {{200, -1.0, 1.0, "x_{vtx}^{rec} (cm)"}});
    registryData.add("mc/reconstruction/collisions/y_vertex", "Reconstructed MC collision vertex;y_{vtx}^{rec} (cm);N_{collisions}", HistType::kTH1F, {{200, -1.0, 1.0, "y_{vtx}^{rec} (cm)"}});
    registryData.add("mc/reconstruction/collisions/xy_vertex", "Reconstructed MC collision vertex;x_{vtx}^{rec} (cm);y_{vtx}^{rec} (cm)", HistType::kTH2F, {{200, -1.0, 1.0, "x_{vtx}^{rec} (cm)"}, {200, -1.0, 1.0, "y_{vtx}^{rec} (cm)"}});
    registryData.add("mc/reconstruction/collisions/n_contributors", "Reconstructed MC collision PV contributors;N_{contributors};N_{collisions}", HistType::kTH1I, {{201, -0.5, 200.5, "N_{contributors}"}});
    registryData.add("mc/reconstruction/collisions/n_selected_tracks", "Selected reconstructed tracks per MC collision;N_{tracks}^{selected};N_{collisions}", HistType::kTH1I, {{201, -0.5, 200.5, "N_{tracks}^{selected}"}});
    registryData.add("mc/reconstruction/collisions/n_mc_matched_tracks", "MC-matched reconstructed tracks per collision;N_{tracks}^{matched};N_{collisions}", HistType::kTH1I, {{201, -0.5, 200.5, "N_{tracks}^{matched}"}});
         
    registryData.add("mc/n_events", "Event counter", HistType::kTH1F, {{1, 0.5, 1.5, "N_{events}"}});

    // Hadrons

    registryData.add("mc/reconstruction/hadrons/tpc_signal_vs_p", "Reconstructed MC tracks;" "#it{p} (GeV/#it{c});" "TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axisMomentumMC, axisTPCSignalMC});
    registryData.add("mc/reconstruction/hadrons/tof_beta_vs_p", "Reconstructed MC tracks with TOF;" "#it{p} (GeV/#it{c});" "#beta_{TOF}", HistType::kTH2F, {axisMomentumMC, axisTOFBetaMC});
    registryData.add("mc/reconstruction/hadrons/rec_hadron_all", "All Tracks (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/hadrons/rec_hadron_tof_match", "All Tracks (Reconstructed) with TOF signal", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/hadrons/mc_rec_hadron_pt", "True Primary Tracks (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/hadrons/mc_sec_hadron_pt", "Secondary Tracks (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/truth/hadrons/mc_gen_hadron_pt", "Transverse-momentum distribution of generated primary charged hadrons;#it{p}_{T} (GeV/#it{c});N_{#pi^{#pm}+K^{#pm}+p+#bar{p}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});

    // Pions

    registryData.add("mc/reconstruction/pions/tpc_signal_vs_p", "True reconstructed #pi^{#pm};" "#it{p} (GeV/#it{c});" "TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axisMomentumMC, axisTPCSignalMC});
    registryData.add("mc/reconstruction/pions/tof_beta_vs_p", "True reconstructed #pi^{#pm} with TOF;" "#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axisMomentumMC, axisTOFBetaMC});
    registryData.add("mc/reconstruction/pions/rec_pion_all", "All Tracks PID'd as Pions", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/mc_rec_pion_pt", "True Primary Pions (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/mc_sec_pion_pt", "Secondary Pions (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/contamination_matrix_pion", "Pion PID Contamination", HistType::kTH2F, {{8000, -4000.5, 3999.5, "PDG Code"}, {120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});

    registryData.add("mc/reconstruction/pions/pos/mc_rec_pion_pos_pt", "Reconstructed Primary #pi^{+} pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/pos/contamination_matrix_pion_pos", "#pi^{+} PID Contamination", HistType::kTH2F, {{8000, -4000.5, 3999.5, "PDG Code"}, {120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/pos/rec_pion_pos_all", "All Tracks PID'd as #pi^{+}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/pos/mc_sec_pion_pos_pt", "Secondary #pi^{+} (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/pos/rec_pion_pos_tof_matched", "Transverse-momentum distribution of TOF-matched reconstructed #pi^{+} in MC;#it{p}_{T} (GeV/#it{c});N_{#pi^{+}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/pos/rec_pion_pos_tpc", "TPC n#sigma_{#pi} vs transverse momentum for reconstructed #pi^{+} tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{#pi}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{#pi}"}});
    registryData.add("mc/reconstruction/pions/pos/rec_pion_pos_tof", "TOF n#sigma_{#pi} vs transverse momentum for reconstructed #pi^{+} tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{#pi}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}^{#pi}"}});

    registryData.add("mc/reconstruction/pions/neg/mc_rec_pion_neg_pt", "Reconstructed Primary #pi^{-} pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/neg/contamination_matrix_pion_neg", "#pi^{-} PID Contamination", HistType::kTH2F, {{8000, -4000.5, 3999.5, "PDG Code"}, {120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/neg/rec_pion_neg_all", "All Tracks PID'd as #pi^{-}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/neg/mc_sec_pion_neg_pt", "Secondary #pi^{-} (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/neg/rec_pion_neg_tof_matched", "Transverse-momentum distribution of TOF-matched reconstructed #pi^{-} in MC;#it{p}_{T} (GeV/#it{c});N_{#pi^{-}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/pions/neg/rec_pion_neg_tpc", "TPC n#sigma_{#pi} vs transverse momentum for reconstructed #pi^{-} tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{#pi}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{#pi}"}});
    registryData.add("mc/reconstruction/pions/neg/rec_pion_neg_tof", "TOF n#sigma_{#pi} vs transverse momentum for reconstructed #pi^{-} tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{#pi}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}^{#pi}"}});

    //Kaons

    registryData.add("mc/reconstruction/kaons/tpc_signal_vs_p", "True reconstructed K^{#pm};" "#it{p} (GeV/#it{c});" "TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axisMomentumMC, axisTPCSignalMC});
    registryData.add("mc/reconstruction/kaons/tof_beta_vs_p", "True reconstructed K^{#pm} with TOF;" "#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axisMomentumMC, axisTOFBetaMC});
    registryData.add("mc/reconstruction/kaons/rec_kaon_all", "All Tracks PID'd as Kaons", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/mc_rec_kaon_pt", "True Primary Kaons (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/mc_sec_kaon_pt", "Secondary Kaons (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/contamination_matrix_kaon", "Kaon PID Contamination", HistType::kTH2F, {{8000, -4000.5, 3999.5, "PDG Code"}, {120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});

    registryData.add("mc/reconstruction/kaons/pos/mc_rec_kaon_pos_pt", "Reconstructed Primary K^{+} pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/pos/contamination_matrix_kaon_pos", "K^{+} PID Contamination", HistType::kTH2F, {{8000, -4000.5, 3999.5, "PDG Code"}, {120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/pos/rec_kaon_pos_all", "All Tracks PID'd as K^{+}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/pos/mc_sec_kaon_pos_pt", "Secondary K^{+} (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/pos/rec_kaon_pos_tof_matched", "Transverse-momentum distribution of TOF-matched reconstructed K^{+} in MC;#it{p}_{T} (GeV/#it{c});N_{K^{+}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/pos/rec_kaon_pos_tpc", "TPC n#sigma_{K} vs transverse momentum for reconstructed K^{+} tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{K}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{K}"}});
    registryData.add("mc/reconstruction/kaons/pos/rec_kaon_pos_tof", "TOF n#sigma_{K} vs transverse momentum for reconstructed K^{+} tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{K}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}^{K}"}});
 
    registryData.add("mc/reconstruction/kaons/neg/rec_kaon_neg_all", "All Tracks PID'd as K^{-}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/neg/mc_sec_kaon_neg_pt", "Secondary K^{-} (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/neg/mc_rec_kaon_neg_pt", "Reconstructed Primary K^{-} pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/neg/contamination_matrix_kaon_neg", "K^{-} PID Contamination", HistType::kTH2F, {{8000, -4000.5, 3999.5, "PDG Code"}, {120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/neg/rec_kaon_neg_tof_matched", "Transverse-momentum distribution of TOF-matched reconstructed K^{-} in MC;#it{p}_{T} (GeV/#it{c});N_{K^{-}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/kaons/neg/rec_kaon_neg_tpc", "TPC n#sigma_{K} vs transverse momentum for reconstructed K^{-} tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{K}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{K}"}});
    registryData.add("mc/reconstruction/kaons/neg/rec_kaon_neg_tof", "TOF n#sigma_{K} vs transverse momentum for reconstructed K^{-} tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{K}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}^{K}"}});

    //Protons

    registryData.add("mc/reconstruction/protons/tpc_signal_vs_p", "True reconstructed p/#bar{p};" "#it{p} (GeV/#it{c});" "TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axisMomentumMC, axisTPCSignalMC});
    registryData.add("mc/reconstruction/protons/tof_beta_vs_p", "True reconstructed p/#bar{p} with TOF;" "#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axisMomentumMC, axisTOFBetaMC});
    registryData.add("mc/reconstruction/protons/rec_proton_all", "All Tracks PID'd as Protons", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/mc_rec_proton_pt", "True Primary Protons (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/mc_sec_proton_pt", "Secondary Protons (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/contamination_matrix_proton", "Proton PID Contamination", HistType::kTH2F, {{8000, -4000.5, 3999.5, "PDG Code"}, {120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});

    registryData.add("mc/reconstruction/protons/pos/mc_rec_proton_pos_pt", "Reconstructed Primary p pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/pos/contamination_matrix_proton_pos", "p PID Contamination", HistType::kTH2F, {{8000, -4000.5, 3999.5, "PDG Code"}, {120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/pos/rec_proton_pos_all", "All Tracks PID'd as p", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/pos/mc_sec_proton_pos_pt", "Secondary p (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/pos/rec_proton_pos_tof_matched", "Transverse-momentum distribution of TOF-matched reconstructed protons in MC;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/pos/rec_proton_pos_tpc", "TPC n#sigma_{p} vs transverse momentum for reconstructed proton tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{p}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{p}"}});
    registryData.add("mc/reconstruction/protons/pos/rec_proton_pos_tof", "TOF n#sigma_{p} vs transverse momentum for reconstructed proton tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{p}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}^{p}"}});

    registryData.add("mc/reconstruction/protons/neg/mc_rec_proton_neg_pt", "Reconstructed Primary #bar{p} pT", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/neg/contamination_matrix_proton_neg", "#bar{p} PID Contamination", HistType::kTH2F, {{8000, -4000.5, 3999.5, "PDG Code"}, {120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/neg/rec_proton_neg_all", "All Tracks PID'd as #bar{p}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/neg/mc_sec_proton_neg_pt", "Secondary #bar{p} (Reconstructed)", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/neg/rec_proton_neg_tof_matched", "Transverse-momentum distribution of TOF-matched reconstructed antiprotons in MC;#it{p}_{T} (GeV/#it{c});N_{#bar{p}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/reconstruction/protons/neg/rec_proton_neg_tpc", "TPC n#sigma_{p} vs transverse momentum for reconstructed antiproton tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{p}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TPC}^{p}"}});
    registryData.add("mc/reconstruction/protons/neg/rec_proton_neg_tof", "TOF n#sigma_{p} vs transverse momentum for reconstructed antiproton tracks in MC;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{p}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -3.0, 3.0, "n#sigma_{TOF}^{p}"}});

    //////////////////////////////////////////////
    //                 MC TRUTH
    //////////////////////////////////////////////

    registryData.add("mc/truth/pions/mc_gen_pion_pt", "Transverse-momentum distribution of generated primary charged pions;#it{p}_{T} (GeV/#it{c});N_{#pi^{+}+#pi^{-}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/truth/pions/pos/mc_gen_pion_pos_pt", "Transverse-momentum distribution of generated primary #pi^{+};#it{p}_{T} (GeV/#it{c});N_{#pi^{+}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/truth/pions/neg/mc_gen_pion_neg_pt", "Transverse-momentum distribution of generated primary #pi^{-};#it{p}_{T} (GeV/#it{c});N_{#pi^{-}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});

    registryData.add("mc/truth/kaons/mc_gen_kaon_pt", "Transverse-momentum distribution of generated primary charged kaons;#it{p}_{T} (GeV/#it{c});N_{K^{+}+K^{-}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/truth/kaons/pos/mc_gen_kaon_pos_pt", "Transverse-momentum distribution of generated primary K^{+};#it{p}_{T} (GeV/#it{c});N_{K^{+}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/truth/kaons/neg/mc_gen_kaon_neg_pt", "Transverse-momentum distribution of generated primary K^{-};#it{p}_{T} (GeV/#it{c});N_{K^{-}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});

    registryData.add("mc/truth/protons/mc_gen_proton_pt", "Transverse-momentum distribution of generated primary protons and antiprotons;#it{p}_{T} (GeV/#it{c});N_{p+#bar{p}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/truth/protons/pos/mc_gen_proton_pos_pt", "Transverse-momentum distribution of generated primary protons;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/truth/protons/neg/mc_gen_proton_neg_pt", "Transverse-momentum distribution of generated primary antiprotons;#it{p}_{T} (GeV/#it{c});N_{#bar{p}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});

    registryData.add("mc/truth/collisions/n_events", "Selected MC-truth collision counter;Selection;N_{collisions}", HistType::kTH1I, {{1, 0.5, 1.5, "Selection"}});
    registryData.add("mc/truth/collisions/z_vertex", "Longitudinal vertex-position distribution of selected MC-truth collisions;z_{vtx}^{MC} (cm);N_{collisions}", HistType::kTH1F, {{200, -20.0, 20.0, "z_{vtx}^{MC} (cm)"}});
    registryData.add("mc/truth/collisions/x_vertex", "Transverse x-vertex distribution of selected MC-truth collisions;x_{vtx}^{MC} (cm);N_{collisions}", HistType::kTH1F, {{200, -1.0, 1.0, "x_{vtx}^{MC} (cm)"}});
    registryData.add("mc/truth/collisions/y_vertex", "Transverse y-vertex distribution of selected MC-truth collisions;y_{vtx}^{MC} (cm);N_{collisions}", HistType::kTH1F, {{200, -1.0, 1.0, "y_{vtx}^{MC} (cm)"}});
    registryData.add("mc/truth/collisions/xy_vertex", "MC-truth collision-vertex distribution in the transverse plane;x_{vtx}^{MC} (cm);y_{vtx}^{MC} (cm)", HistType::kTH2F, {{200, -1.0, 1.0, "x_{vtx}^{MC} (cm)"}, {200, -1.0, 1.0, "y_{vtx}^{MC} (cm)"}});
    registryData.add("mc/truth/collisions/n_selected_hadrons", "Multiplicity of selected primary charged hadrons per MC-truth collision;N_{#pi^{#pm}+K^{#pm}+p+#bar{p}};N_{collisions}", HistType::kTH1I, {axisTruthMultiplicity});

    registryData.add("mc/truth/hadrons/mc_gen_hadron_eta", "Pseudorapidity distribution of generated primary #pi^{#pm}, K^{#pm}, p and #bar{p};#eta;N_{particles}", HistType::kTH1F, {axisTruthEta});
    registryData.add("mc/truth/hadrons/mc_gen_hadron_phi", "Azimuthal-angle distribution of generated primary #pi^{#pm}, K^{#pm}, p and #bar{p};#varphi;N_{particles}", HistType::kTH1F, {axisTruthPhi});
    registryData.add("mc/truth/hadrons/mc_gen_hadron_eta_phi", "Pseudorapidity vs azimuthal angle for generated primary #pi^{#pm}, K^{#pm}, p and #bar{p};#varphi;#eta", HistType::kTH2F, {axisTruthPhi, axisTruthEta});
    registryData.add("mc/truth/hadrons/mc_gen_hadron_phi_vs_pt", "Azimuthal angle vs transverse momentum for generated primary #pi^{#pm}, K^{#pm}, p and #bar{p};#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axisTruthPt, axisTruthPhi});

    // Pions

    registryData.add("mc/truth/pions/mc_gen_pion_eta", "Pseudorapidity distribution of generated primary charged pions;#eta;N_{#pi^{+}+#pi^{-}}", HistType::kTH1F, {axisTruthEta});
    registryData.add("mc/truth/pions/mc_gen_pion_phi", "Azimuthal-angle distribution of generated primary charged pions;#varphi;N_{#pi^{+}+#pi^{-}}", HistType::kTH1F, {axisTruthPhi});
    registryData.add("mc/truth/pions/mc_gen_pion_eta_phi", "Pseudorapidity vs azimuthal angle for generated primary charged pions;#varphi;#eta", HistType::kTH2F, {axisTruthPhi, axisTruthEta});
    registryData.add("mc/truth/pions/mc_gen_pion_phi_vs_pt", "Azimuthal angle vs transverse momentum for generated primary charged pions;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axisTruthPt, axisTruthPhi});

    registryData.add("mc/truth/pions/pos/mc_gen_pion_pos_eta", "Pseudorapidity distribution of generated primary #pi^{+};#eta;N_{#pi^{+}}", HistType::kTH1F, {axisTruthEta});
    registryData.add("mc/truth/pions/pos/mc_gen_pion_pos_phi", "Azimuthal-angle distribution of generated primary #pi^{+};#varphi;N_{#pi^{+}}", HistType::kTH1F, {axisTruthPhi});
    registryData.add("mc/truth/pions/pos/mc_gen_pion_pos_eta_phi", "Pseudorapidity vs azimuthal angle for generated primary #pi^{+};#varphi;#eta", HistType::kTH2F, {axisTruthPhi, axisTruthEta});
    registryData.add("mc/truth/pions/pos/mc_gen_pion_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for generated primary #pi^{+};#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axisTruthPt, axisTruthPhi});

    registryData.add("mc/truth/pions/neg/mc_gen_pion_neg_eta", "Pseudorapidity distribution of generated primary #pi^{-};#eta;N_{#pi^{-}}", HistType::kTH1F, {axisTruthEta});
    registryData.add("mc/truth/pions/neg/mc_gen_pion_neg_phi", "Azimuthal-angle distribution of generated primary #pi^{-};#varphi;N_{#pi^{-}}", HistType::kTH1F, {axisTruthPhi});
    registryData.add("mc/truth/pions/neg/mc_gen_pion_neg_eta_phi", "Pseudorapidity vs azimuthal angle for generated primary #pi^{-};#varphi;#eta", HistType::kTH2F, {axisTruthPhi, axisTruthEta});
    registryData.add("mc/truth/pions/neg/mc_gen_pion_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for generated primary #pi^{-};#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axisTruthPt, axisTruthPhi});

    // Kaons

    registryData.add("mc/truth/kaons/mc_gen_kaon_eta", "Pseudorapidity distribution of generated primary charged kaons;#eta;N_{K^{+}+K^{-}}", HistType::kTH1F, {axisTruthEta});
    registryData.add("mc/truth/kaons/mc_gen_kaon_phi", "Azimuthal-angle distribution of generated primary charged kaons;#varphi;N_{K^{+}+K^{-}}", HistType::kTH1F, {axisTruthPhi});
    registryData.add("mc/truth/kaons/mc_gen_kaon_eta_phi", "Pseudorapidity vs azimuthal angle for generated primary charged kaons;#varphi;#eta", HistType::kTH2F, {axisTruthPhi, axisTruthEta});
    registryData.add("mc/truth/kaons/mc_gen_kaon_phi_vs_pt", "Azimuthal angle vs transverse momentum for generated primary charged kaons;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axisTruthPt, axisTruthPhi});

    registryData.add("mc/truth/kaons/pos/mc_gen_kaon_pos_eta", "Pseudorapidity distribution of generated primary K^{+};#eta;N_{K^{+}}", HistType::kTH1F, {axisTruthEta});
    registryData.add("mc/truth/kaons/pos/mc_gen_kaon_pos_phi", "Azimuthal-angle distribution of generated primary K^{+};#varphi;N_{K^{+}}", HistType::kTH1F, {axisTruthPhi});
    registryData.add("mc/truth/kaons/pos/mc_gen_kaon_pos_eta_phi", "Pseudorapidity vs azimuthal angle for generated primary K^{+};#varphi;#eta", HistType::kTH2F, {axisTruthPhi, axisTruthEta});
    registryData.add("mc/truth/kaons/pos/mc_gen_kaon_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for generated primary K^{+};#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axisTruthPt, axisTruthPhi});

    registryData.add("mc/truth/kaons/neg/mc_gen_kaon_neg_eta", "Pseudorapidity distribution of generated primary K^{-};#eta;N_{K^{-}}", HistType::kTH1F, {axisTruthEta});
    registryData.add("mc/truth/kaons/neg/mc_gen_kaon_neg_phi", "Azimuthal-angle distribution of generated primary K^{-};#varphi;N_{K^{-}}", HistType::kTH1F, {axisTruthPhi});
    registryData.add("mc/truth/kaons/neg/mc_gen_kaon_neg_eta_phi", "Pseudorapidity vs azimuthal angle for generated primary K^{-};#varphi;#eta", HistType::kTH2F, {axisTruthPhi, axisTruthEta});
    registryData.add("mc/truth/kaons/neg/mc_gen_kaon_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for generated primary K^{-};#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axisTruthPt, axisTruthPhi});

    // Protons and antiprotons

    registryData.add("mc/truth/protons/mc_gen_proton_eta", "Pseudorapidity distribution of generated primary protons and antiprotons;#eta;N_{p+#bar{p}}", HistType::kTH1F, {axisTruthEta});
    registryData.add("mc/truth/protons/mc_gen_proton_phi", "Azimuthal-angle distribution of generated primary protons and antiprotons;#varphi;N_{p+#bar{p}}", HistType::kTH1F, {axisTruthPhi});
    registryData.add("mc/truth/protons/mc_gen_proton_eta_phi", "Pseudorapidity vs azimuthal angle for generated primary protons and antiprotons;#varphi;#eta", HistType::kTH2F, {axisTruthPhi, axisTruthEta});
    registryData.add("mc/truth/protons/mc_gen_proton_phi_vs_pt", "Azimuthal angle vs transverse momentum for generated primary protons and antiprotons;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axisTruthPt, axisTruthPhi});

    registryData.add("mc/truth/protons/pos/mc_gen_proton_pos_eta", "Pseudorapidity distribution of generated primary protons;#eta;N_{p}", HistType::kTH1F, {axisTruthEta});
    registryData.add("mc/truth/protons/pos/mc_gen_proton_pos_phi", "Azimuthal-angle distribution of generated primary protons;#varphi;N_{p}", HistType::kTH1F, {axisTruthPhi});
    registryData.add("mc/truth/protons/pos/mc_gen_proton_pos_eta_phi", "Pseudorapidity vs azimuthal angle for generated primary protons;#varphi;#eta", HistType::kTH2F, {axisTruthPhi, axisTruthEta});
    registryData.add("mc/truth/protons/pos/mc_gen_proton_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for generated primary protons;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axisTruthPt, axisTruthPhi});

    registryData.add("mc/truth/protons/neg/mc_gen_proton_neg_eta", "Pseudorapidity distribution of generated primary antiprotons;#eta;N_{#bar{p}}", HistType::kTH1F, {axisTruthEta});
    registryData.add("mc/truth/protons/neg/mc_gen_proton_neg_phi", "Azimuthal-angle distribution of generated primary antiprotons;#varphi;N_{#bar{p}}", HistType::kTH1F, {axisTruthPhi});
    registryData.add("mc/truth/protons/neg/mc_gen_proton_neg_eta_phi", "Pseudorapidity vs azimuthal angle for generated primary antiprotons;#varphi;#eta", HistType::kTH2F, {axisTruthPhi, axisTruthEta});
    registryData.add("mc/truth/protons/neg/mc_gen_proton_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for generated primary antiprotons;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axisTruthPt, axisTruthPhi});
        
    //////////////////////////////////////////////
    //         MC DETECTOR-LEVEL JETS
    //////////////////////////////////////////////

    registryData.add("mc/jets/detector/n_events", "Selected events with detector-level charged jets;Event counter;N_{events}", HistType::kTH1I, {{1, 0.5, 1.5, "Event counter"}});
    registryData.add("mc/jets/detector/jet_pt", "Transverse-momentum distribution of detector-level charged jets;#it{p}_{T,jet}^{det} (GeV/#it{c});N_{jets}", HistType::kTH1F, {{200, 0.0, 100.0, "#it{p}_{T,jet}^{det} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/jet_eta", "Pseudorapidity distribution of detector-level charged jets;#eta_{jet}^{det};N_{jets}", HistType::kTH1F, {{100, -1.0, 1.0, "#eta_{jet}^{det}"}});
    registryData.add("mc/jets/detector/jet_phi", "Azimuthal-angle distribution of detector-level charged jets;#varphi_{jet}^{det};N_{jets}", HistType::kTH1F, {{72, 0.0, TwoPI, "#varphi_{jet}^{det}"}});
    registryData.add("mc/jets/detector/jet_n_constituents", "Constituent multiplicity of detector-level charged jets;N_{constituents}^{det};N_{jets}", HistType::kTH1I, {{101, -0.5, 100.5, "N_{constituents}^{det}"}});
    registryData.add("mc/jets/detector/jet_area", "Area distribution of detector-level charged jets;A_{jet}^{det};N_{jets}", HistType::kTH1F, {{100, 0.0, 1.5, "A_{jet}^{det}"}});
    registryData.add("mc/jets/detector/jet_n_selected_constituents", "Selected-track multiplicity in detector-level charged jets;N_{constituents}^{selected};N_{jets}", HistType::kTH1I, {{101, -0.5, 100.5, "N_{constituents}^{selected}"}});

    // Pions

    registryData.add("mc/jets/detector/pions/pos/pion_jet_pos_tpc", "TPC n#sigma_{#pi} vs transverse momentum for #pi^{+} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{#pi}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TPC}^{#pi}"}});
    registryData.add("mc/jets/detector/pions/pos/pion_jet_pos_tof", "TOF n#sigma_{#pi} vs transverse momentum for #pi^{+} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{#pi}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TOF}^{#pi}"}});
    registryData.add("mc/jets/detector/pions/pos/pion_jet_pos_pt", "Transverse-momentum distribution of #pi^{+} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});N_{#pi^{+}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/pions/pos/pion_jet_pos_eta", "Pseudorapidity distribution of #pi^{+} candidates in detector-level MC jets;#eta;N_{#pi^{+}}", HistType::kTH1F, {{80, -0.8, 0.8, "#eta"}});
    registryData.add("mc/jets/detector/pions/pos/pion_jet_pos_dcaxy", "DCA_{xy} distribution of #pi^{+} candidates in detector-level MC jets;DCA_{xy} (cm);N_{#pi^{+}}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{xy} (cm)"}});
    registryData.add("mc/jets/detector/pions/pos/pion_jet_pos_dcaz", "DCA_{z} distribution of #pi^{+} candidates in detector-level MC jets;DCA_{z} (cm);N_{#pi^{+}}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{z} (cm)"}});
    registryData.add("mc/jets/detector/pions/pos/pion_jet_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for #pi^{+} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});

    registryData.add("mc/jets/detector/pions/neg/pion_jet_neg_tpc", "TPC n#sigma_{#pi} vs transverse momentum for #pi^{-} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{#pi}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TPC}^{#pi}"}});
    registryData.add("mc/jets/detector/pions/neg/pion_jet_neg_tof", "TOF n#sigma_{#pi} vs transverse momentum for #pi^{-} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{#pi}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TOF}^{#pi}"}});
    registryData.add("mc/jets/detector/pions/neg/pion_jet_neg_pt", "Transverse-momentum distribution of #pi^{-} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});N_{#pi^{-}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/pions/neg/pion_jet_neg_eta", "Pseudorapidity distribution of #pi^{-} candidates in detector-level MC jets;#eta;N_{#pi^{-}}", HistType::kTH1F, {{80, -0.8, 0.8, "#eta"}});
    registryData.add("mc/jets/detector/pions/neg/pion_jet_neg_dcaxy", "DCA_{xy} distribution of #pi^{-} candidates in detector-level MC jets;DCA_{xy} (cm);N_{#pi^{-}}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{xy} (cm)"}});
    registryData.add("mc/jets/detector/pions/neg/pion_jet_neg_dcaz", "DCA_{z} distribution of #pi^{-} candidates in detector-level MC jets;DCA_{z} (cm);N_{#pi^{-}}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{z} (cm)"}});
    registryData.add("mc/jets/detector/pions/neg/pion_jet_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for #pi^{-} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});

    // Kaons

    registryData.add("mc/jets/detector/kaons/pos/kaon_jet_pos_tpc", "TPC n#sigma_{K} vs transverse momentum for K^{+} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{K}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TPC}^{K}"}});
    registryData.add("mc/jets/detector/kaons/pos/kaon_jet_pos_tof", "TOF n#sigma_{K} vs transverse momentum for K^{+} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{K}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TOF}^{K}"}});
    registryData.add("mc/jets/detector/kaons/pos/kaon_jet_pos_pt", "Transverse-momentum distribution of K^{+} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});N_{K^{+}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/kaons/pos/kaon_jet_pos_eta", "Pseudorapidity distribution of K^{+} candidates in detector-level MC jets;#eta;N_{K^{+}}", HistType::kTH1F, {{80, -0.8, 0.8, "#eta"}});
    registryData.add("mc/jets/detector/kaons/pos/kaon_jet_pos_dcaxy", "DCA_{xy} distribution of K^{+} candidates in detector-level MC jets;DCA_{xy} (cm);N_{K^{+}}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{xy} (cm)"}});
    registryData.add("mc/jets/detector/kaons/pos/kaon_jet_pos_dcaz", "DCA_{z} distribution of K^{+} candidates in detector-level MC jets;DCA_{z} (cm);N_{K^{+}}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{z} (cm)"}});
    registryData.add("mc/jets/detector/kaons/pos/kaon_jet_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for K^{+} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});

    registryData.add("mc/jets/detector/kaons/neg/kaon_jet_neg_tpc", "TPC n#sigma_{K} vs transverse momentum for K^{-} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{K}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TPC}^{K}"}});
    registryData.add("mc/jets/detector/kaons/neg/kaon_jet_neg_tof", "TOF n#sigma_{K} vs transverse momentum for K^{-} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{K}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TOF}^{K}"}});
    registryData.add("mc/jets/detector/kaons/neg/kaon_jet_neg_pt", "Transverse-momentum distribution of K^{-} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});N_{K^{-}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/kaons/neg/kaon_jet_neg_eta", "Pseudorapidity distribution of K^{-} candidates in detector-level MC jets;#eta;N_{K^{-}}", HistType::kTH1F, {{80, -0.8, 0.8, "#eta"}});
    registryData.add("mc/jets/detector/kaons/neg/kaon_jet_neg_dcaxy", "DCA_{xy} distribution of K^{-} candidates in detector-level MC jets;DCA_{xy} (cm);N_{K^{-}}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{xy} (cm)"}});
    registryData.add("mc/jets/detector/kaons/neg/kaon_jet_neg_dcaz", "DCA_{z} distribution of K^{-} candidates in detector-level MC jets;DCA_{z} (cm);N_{K^{-}}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{z} (cm)"}});
    registryData.add("mc/jets/detector/kaons/neg/kaon_jet_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for K^{-} candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});

    // Protons and antiprotons

    registryData.add("mc/jets/detector/protons/pos/proton_jet_pos_tpc", "TPC n#sigma_{p} vs transverse momentum for proton candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{p}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TPC}^{p}"}});
    registryData.add("mc/jets/detector/protons/pos/proton_jet_pos_tof", "TOF n#sigma_{p} vs transverse momentum for proton candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{p}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TOF}^{p}"}});
    registryData.add("mc/jets/detector/protons/pos/proton_jet_pos_pt", "Transverse-momentum distribution of proton candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/protons/pos/proton_jet_pos_eta", "Pseudorapidity distribution of proton candidates in detector-level MC jets;#eta;N_{p}", HistType::kTH1F, {{80, -0.8, 0.8, "#eta"}});
    registryData.add("mc/jets/detector/protons/pos/proton_jet_pos_dcaxy", "DCA_{xy} distribution of proton candidates in detector-level MC jets;DCA_{xy} (cm);N_{p}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{xy} (cm)"}});
    registryData.add("mc/jets/detector/protons/pos/proton_jet_pos_dcaz", "DCA_{z} distribution of proton candidates in detector-level MC jets;DCA_{z} (cm);N_{p}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{z} (cm)"}});
    registryData.add("mc/jets/detector/protons/pos/proton_jet_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for proton candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});

    registryData.add("mc/jets/detector/protons/neg/antiproton_jet_neg_tpc", "TPC n#sigma_{p} vs transverse momentum for antiproton candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{p}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TPC}^{p}"}});
    registryData.add("mc/jets/detector/protons/neg/antiproton_jet_neg_tof", "TOF n#sigma_{p} vs transverse momentum for antiproton candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{p}", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -5.0, 5.0, "n#sigma_{TOF}^{p}"}});
    registryData.add("mc/jets/detector/protons/neg/antiproton_jet_neg_pt", "Transverse-momentum distribution of antiproton candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});N_{#bar{p}}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/protons/neg/antiproton_jet_neg_eta", "Pseudorapidity distribution of antiproton candidates in detector-level MC jets;#eta;N_{#bar{p}}", HistType::kTH1F, {{80, -0.8, 0.8, "#eta"}});
    registryData.add("mc/jets/detector/protons/neg/antiproton_jet_neg_dcaxy", "DCA_{xy} distribution of antiproton candidates in detector-level MC jets;DCA_{xy} (cm);N_{#bar{p}}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{xy} (cm)"}});
    registryData.add("mc/jets/detector/protons/neg/antiproton_jet_neg_dcaz", "DCA_{z} distribution of antiproton candidates in detector-level MC jets;DCA_{z} (cm);N_{#bar{p}}", HistType::kTH1F, {{200, -0.2, 0.2, "DCA_{z} (cm)"}});
    registryData.add("mc/jets/detector/protons/neg/antiproton_jet_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for antiproton candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}, {72, 0.0, TwoPI, "#varphi"}});

    //////////////////////////////////////////////
    //              MCD JET CONE
    //////////////////////////////////////////////

    registryData.add("mc/jets/detector/cone/all_tracks_2d", "Selected tracks in MCD jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone});
    registryData.add("mc/jets/detector/cone/pions_2d", "Reconstructed #pi in MCD jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone});
    registryData.add("mc/jets/detector/cone/kaons_2d", "Reconstructed K in MCD jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone});
    registryData.add("mc/jets/detector/cone/protons_2d", "Reconstructed p/#bar{p} in MCD jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone});

    registryData.add("mc/jets/detector/cone/all_tracks_3d", "Selected tracks in MCD jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone, axisMCDPtTrackCone});
    registryData.add("mc/jets/detector/cone/pions_3d", "Reconstructed #pi in MCD jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone, axisMCDPtTrackCone});
    registryData.add("mc/jets/detector/cone/kaons_3d", "Reconstructed K in MCD jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone, axisMCDPtTrackCone});
    registryData.add("mc/jets/detector/cone/protons_3d", "Reconstructed p/#bar{p} in MCD jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone, axisMCDPtTrackCone});

    //////////////////////////////////////////////
    //           MCD JET COMPOSITION
    //////////////////////////////////////////////

    registryData.add("mc/jets/detector/composition/n_pions_per_jet", "Multiplicity of PID-selected pion candidates in detector-level MC jets;N_{#pi^{+}+#pi^{-}};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/composition/n_kaons_per_jet", "Multiplicity of PID-selected kaon candidates in detector-level MC jets;N_{K^{+}+K^{-}};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/composition/n_protons_per_jet", "Multiplicity of PID-selected proton and antiproton candidates in detector-level MC jets;N_{p+#bar{p}};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});

    registryData.add("mc/jets/detector/composition/n_pions_vs_n_kaons", "Kaon vs pion multiplicity in detector-level MC jets;N_{#pi^{+}+#pi^{-}};N_{K^{+}+K^{-}};N_{jets}", HistType::kTH2I, {axisMCDNPions, axisMCDNKaons});
    registryData.add("mc/jets/detector/composition/n_pions_vs_n_protons", "Proton vs pion multiplicity in detector-level MC jets;N_{#pi^{+}+#pi^{-}};N_{p+#bar{p}};N_{jets}", HistType::kTH2I, {axisMCDNPions, axisMCDNProtons});
    registryData.add("mc/jets/detector/composition/n_kaons_vs_n_protons", "Proton vs kaon multiplicity in detector-level MC jets;N_{K^{+}+K^{-}};N_{p+#bar{p}};N_{jets}", HistType::kTH2I, {axisMCDNKaons, axisMCDNProtons});
    registryData.add("mc/jets/detector/composition/n_pions_n_kaons_n_protons", "Joint identified-particle multiplicity in detector-level MC jets;N_{#pi^{+}+#pi^{-}};N_{K^{+}+K^{-}};N_{p+#bar{p}}", HistType::kTH3F, {axisMCDNPions, axisMCDNKaons, axisMCDNProtons});

    registryData.add("mc/jets/detector/composition/pions/pos/n_pions_per_jet", "Multiplicity of PID-selected #pi^{+} candidates in detector-level MC jets;N_{#pi^{+}};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/composition/pions/neg/n_pions_per_jet", "Multiplicity of PID-selected #pi^{-} candidates in detector-level MC jets;N_{#pi^{-}};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});

    registryData.add("mc/jets/detector/composition/kaons/pos/n_kaons_per_jet", "Multiplicity of PID-selected K^{+} candidates in detector-level MC jets;N_{K^{+}};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/composition/kaons/neg/n_kaons_per_jet", "Multiplicity of PID-selected K^{-} candidates in detector-level MC jets;N_{K^{-}};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/composition/protons/pos/n_protons_per_jet", "Multiplicity of PID-selected proton candidates in detector-level MC jets;N_{p};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/composition/protons/neg/n_antiprotons_per_jet", "Multiplicity of PID-selected antiproton candidates in detector-level MC jets;N_{#bar{p}};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});

    registryData.add("mc/jets/detector/pions/pion_pt", "Transverse-momentum distribution of PID-selected pion candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});N_{#pi candidates}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/kaons/kaon_pt", "Transverse-momentum distribution of PID-selected kaon candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});N_{K candidates}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/protons/proton_pt", "Transverse-momentum distribution of PID-selected proton and antiproton candidates in detector-level MC jets;#it{p}_{T} (GeV/#it{c});N_{p+#bar{p} candidates}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
        
    //Flavor
    registryData.add("mc/jets/detector/flavor/quark_jet_pt", "Quark Jets;#it{p}_{T,jet}^{det} (GeV/#it{c});N_{jets}", HistType::kTH1F, {{200, 0.0, 100.0, "#it{p}_{T,jet}^{det} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/flavor/gluon_jet_pt", "Gluon Jets;#it{p}_{T,jet}^{det} (GeV/#it{c});N_{jets}", HistType::kTH1F, {{200, 0.0, 100.0, "#it{p}_{T,jet}^{det} (GeV/#it{c})"}});

    registryData.add("mc/jets/detector/flavor/quark_pions_per_jet", "Reconstructed #pi in Quark jets;N_{#pi};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/flavor/gluon_pions_per_jet", "Reconstructed #pi in Gluon jets;N_{#pi};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    
    registryData.add("mc/jets/detector/flavor/quark_kaons_per_jet", "Reconstructed K in Quark jets;N_{K};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/flavor/gluon_kaons_per_jet", "Reconstructed K in Gluon jets;N_{K};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    
    registryData.add("mc/jets/detector/flavor/quark_protons_per_jet", "Reconstructed p/#bar{p} in Quark jets;N_{p+#bar{p}};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/flavor/gluon_protons_per_jet", "Reconstructed p/#bar{p} in Gluon jets;N_{p+#bar{p}};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});

    //////////////////////////////////////////////
    //              MCD UE
    //////////////////////////////////////////////

    registryData.add("mc/jets/detector/ue/all_tracks_2d", "Selected tracks in MCD perpendicular cones;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone});
    registryData.add("mc/jets/detector/ue/pions_2d", "Reconstructed #pi in MCD perpendicular cones;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone});
    registryData.add("mc/jets/detector/ue/kaons_2d", "Reconstructed K in MCD perpendicular cones;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone});
    registryData.add("mc/jets/detector/ue/protons_2d", "Reconstructed p/#bar{p} in MCD perpendicular cones;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisMCDDeltaEtaCone, axisMCDDeltaPhiCone});

    registryData.add("mc/jets/detector/ue/n_tracks_per_jet", "Selected tracks in both MCD perpendicular cones;N_{tracks}^{UE1+UE2};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/ue/n_pions_per_jet", "Reconstructed #pi in both MCD perpendicular cones;N_{#pi}^{UE1+UE2};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/ue/n_kaons_per_jet", "Reconstructed K in both MCD perpendicular cones;N_{K}^{UE1+UE2};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});
    registryData.add("mc/jets/detector/ue/n_protons_per_jet", "Reconstructed p/#bar{p} in both MCD perpendicular cones;N_{p+#bar{p}}^{UE1+UE2};N_{jets}", HistType::kTH1I, {axisMCDMultiplicity});

    registryData.add("mc/jets/detector/ue/pion_pt", "Reconstructed #pi in MCD perpendicular cones;#it{p}_{T} (GeV/#it{c});N_{tracks}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/ue/kaon_pt", "Reconstructed K in MCD perpendicular cones;#it{p}_{T} (GeV/#it{c});N_{tracks}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryData.add("mc/jets/detector/ue/proton_pt", "Reconstructed p/#bar{p} in MCD perpendicular cones;#it{p}_{T} (GeV/#it{c});N_{tracks}", HistType::kTH1F, {{120, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"}});  
  }

  void getPerpendicularDirections(const TVector3& p, TVector3& u1, TVector3& u2)
  {
    double const threshold = 1e-9;

    if (p.Pt() < threshold) {
      u1.SetXYZ(0, 0, 0);
      u2.SetXYZ(0, 0, 0);
      return;
    }

    double pt = p.Pt();
    double eta = p.Eta();
    double phi = p.Phi();

    u1.SetPtEtaPhi(pt, eta, phi + PIHalf);
    u2.SetPtEtaPhi(pt, eta, phi - PIHalf);
  }

  template <typename TrackIts>
  bool hasITSLayerHit(const TrackIts& track, int layer)
  {
    return TESTBIT(track.itsClusterMap(), layer - 1);
  }

  struct PidResult {
    bool isPion, isKaon, isProton;
  };

  template <typename TrackType>
  PidResult getPid(const TrackType& track)
  {
    constexpr int ClosestMatch = 0;
    constexpr int ExclusiveMatch = 1;
    constexpr int RejectionBased = 2;

    double const buffer = 999.0;
    double pt = track.pt();

    double dPi = 0;
    double dKa = 0;
    double dPr = 0;

    if (pt < cfg.ptThresholdPion) {
      dPi = std::abs(track.tpcNSigmaPi());
    } else {
      if (track.hasTOF()) {
        dPi = std::hypot(track.tofNSigmaPi(), track.tpcNSigmaPi());
      } else {
        dPi = buffer;
      }
    }

    if (pt < cfg.ptThresholdKaon) {
      dKa = std::abs(track.tpcNSigmaKa());
    } else {
      if (track.hasTOF()) {
        dKa = std::hypot(track.tofNSigmaKa(), track.tpcNSigmaKa());
      } else {
        dKa = buffer;
      }
    }

    if (pt < cfg.ptThresholdProton) {
      dPr = std::abs(track.tpcNSigmaPr());
    } else {
      if (track.hasTOF()) {
        dPr = std::hypot(track.tofNSigmaPr(), track.tpcNSigmaPr());
      } else {
        dPr = buffer;
      }
    }

    bool isPiMatch = (dPi <= cfg.nSigmaCut);
    bool isKaMatch = (dKa <= cfg.nSigmaCut);
    bool isPrMatch = (dPr <= cfg.nSigmaCut);

    PidResult res{.isPion = false, .isKaon = false, .isProton = false};

    if (cfg.pidMethod == ClosestMatch) {
      if (isPiMatch && dPi < dKa && dPi < dPr) {
        res.isPion = true;
      } else if (isKaMatch && dKa < dPi && dKa < dPr) {
        res.isKaon = true;
      } else if (isPrMatch && dPr < dPi && dPr < dKa) {
        res.isProton = true;
      }
    } else if (cfg.pidMethod == ExclusiveMatch) {
      if (isPiMatch && !isKaMatch && !isPrMatch) {
        res.isPion = true;
      } else if (isKaMatch && !isPiMatch && !isPrMatch) {
        res.isKaon = true;
      } else if (isPrMatch && !isPiMatch && !isKaMatch) {
        res.isProton = true;
      }
    } else if (cfg.pidMethod == RejectionBased) {
      if (isPiMatch && dKa > cfg.rejectionSigma && dPr > cfg.rejectionSigma) {
        res.isPion = true;
      } else if (isKaMatch && dPi > cfg.rejectionSigma && dPr > cfg.rejectionSigma) {
        res.isKaon = true;
      } else if (isPrMatch && dPi > cfg.rejectionSigma && dKa > cfg.rejectionSigma) {
        res.isProton = true;
      }
    }

    if (res.isPion && (pt < cfg.minPtPion || pt > cfg.maxPtPion)) {
      res.isPion = false;
    }
    if (res.isKaon && (pt < cfg.minPtKaon || pt > cfg.maxPtKaon)) {
      res.isKaon = false;
    }
    if (res.isProton && (pt < cfg.minPtProton || pt > cfg.maxPtProton)) {
      res.isProton = false;
    }

    return res;
  }

  template <typename TrackType>
  bool passedTrackSelection(const TrackType& track)
  {
    if (requirePvContributor && !(track.isPVContributor())) {
      return false;
    }
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }
    if ((!hasITSLayerHit(track, 1)) && (!hasITSLayerHit(track, 2)) && (!hasITSLayerHit(track, 3))) {
      return false;
    }
    if (track.itsNCls() < minItsNclusters) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < minTpcNcrossedRows) {
      return false;
    }
    if (track.tpcChi2NCl() < minChiSquareTpc || track.tpcChi2NCl() > maxChiSquareTpc) {
      return false;
    }
    if (track.itsChi2NCl() > maxChiSquareIts) {
      return false;
    }
    if (track.eta() < minEta || track.eta() > maxEta) {
      return false;
    }
    if (track.pt() < minPt || track.pt() > maxPt) {
      return false;
    }
    if (std::abs(track.dcaXY()) > maxDcaxy || std::abs(track.dcaZ()) > maxDcaz) {
      return false;
    }
    return true;
  }

  void processPureTracks(StandardEvents::iterator const& collision, HadronTracks const& globalTracks)
  {
    if (!collision.sel8() || std::abs(collision.posZ()) > zVtx) {
      return;
    }
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return;
    }

    int nSelectedTracks = 0;

    registryData.fill(HIST("data/pure/collisions/z_vertex"), collision.posZ());
    registryData.fill(HIST("data/pure/collisions/x_vertex"), collision.posX());
    registryData.fill(HIST("data/pure/collisions/y_vertex"), collision.posY());
    registryData.fill(HIST("data/pure/collisions/xy_vertex"), collision.posX(), collision.posY());
    registryData.fill(HIST("data/pure/collisions/n_contributors"), collision.numContrib());

    for (auto const& track : globalTracks) {
      if (!passedTrackSelection(track)) {
        continue;
      }

      ++nSelectedTracks;

      const double p = track.p();
      registryData.fill(HIST("data/pure/tpc_signal_vs_p"), p, track.tpcSignal());
      if (track.hasTOF()) {
        registryData.fill(HIST("data/pure/tof_beta_vs_p"), p, track.beta());
      }

      double pt = track.pt();
      double eta = track.eta();
      double dcaxy = track.dcaXY();
      double dcaz = track.dcaZ();
      int charge = track.sign();
      double phi = track.phi();

      PidResult pid = getPid(track);

      if (pid.isPion) {
        registryData.fill(HIST("data/pure/pions/pion_tpc_signal_vs_p"), p, track.tpcSignal());
        registryData.fill(HIST("data/pure/pions/pion_pure_tpc"), pt, track.tpcNSigmaPi());
        registryData.fill(HIST("data/pure/pions/pion_pure_eta_phi"), phi, eta);
        if (track.hasTOF()) {
          registryData.fill(HIST("data/pure/pions/pion_tof_beta_vs_p"), p, track.beta());
          registryData.fill(HIST("data/pure/pions/pion_pure_tof"), pt, track.tofNSigmaPi());
        }
        registryData.fill(HIST("data/pure/pions/pion_pure_pt"), pt);
        registryData.fill(HIST("data/pure/pions/pion_pure_eta"), eta);
        registryData.fill(HIST("data/pure/pions/pion_pure_dcaxy"), dcaxy);
        registryData.fill(HIST("data/pure/pions/pion_pure_dcaz"), dcaz);

        if (charge > 0) {
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_tpc"), pt, track.tpcNSigmaPi());
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_tof"), pt, track.tofNSigmaPi());
          }
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_pt"), pt);
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_eta"), eta);
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_dcaxy"), dcaxy);
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_dcaz"), dcaz);
        } else {
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_tpc"), pt, track.tpcNSigmaPi());
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_tof"), pt, track.tofNSigmaPi());
          }
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_pt"), pt);
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_eta"), eta);
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_dcaxy"), dcaxy);
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_dcaz"), dcaz);
        }
      }

      if (pid.isKaon) {
        registryData.fill(HIST("data/pure/kaons/kaon_tpc_signal_vs_p"), track.p(), track.tpcSignal());
        registryData.fill(HIST("data/pure/kaons/kaon_pure_tpc"), pt, track.tpcNSigmaKa());
        registryData.fill(HIST("data/pure/kaons/kaon_pure_eta_phi"), phi, eta);
        if (track.hasTOF()) {
          registryData.fill(HIST("data/pure/kaons/kaon_tof_beta_vs_p"), p, track.beta());
          registryData.fill(HIST("data/pure/kaons/kaon_pure_tof"), pt, track.tofNSigmaKa());
        }
        registryData.fill(HIST("data/pure/kaons/kaon_pure_pt"), pt);
        registryData.fill(HIST("data/pure/kaons/kaon_pure_eta"), eta);
        registryData.fill(HIST("data/pure/kaons/kaon_pure_dcaxy"), dcaxy);
        registryData.fill(HIST("data/pure/kaons/kaon_pure_dcaz"), dcaz);

        if (charge > 0) {
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_tpc"), pt, track.tpcNSigmaKa());
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_tof"), pt, track.tofNSigmaKa());
          }
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_pt"), pt);
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_eta"), eta);
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_dcaxy"), dcaxy);
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_dcaz"), dcaz);
        } else {
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_tpc"), pt, track.tpcNSigmaKa());
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_tof"), pt, track.tofNSigmaKa());
          }
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_pt"), pt);
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_eta"), eta);
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_dcaxy"), dcaxy);
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_dcaz"), dcaz);
        }
      }

      if (pid.isProton) {
        registryData.fill(HIST("data/pure/protons/proton_tpc_signal_vs_p"), track.p(), track.tpcSignal());
        registryData.fill(HIST("data/pure/protons/proton_pure_tpc"), pt, track.tpcNSigmaPr());
        registryData.fill(HIST("data/pure/protons/proton_pure_eta_phi"), phi, eta);
        if (track.hasTOF()) {
          registryData.fill(HIST("data/pure/protons/proton_tof_beta_vs_p"), p, track.beta());
          registryData.fill(HIST("data/pure/protons/proton_pure_tof"), pt, track.tofNSigmaPr());
        }
        registryData.fill(HIST("data/pure/protons/proton_pure_pt"), pt);
        registryData.fill(HIST("data/pure/protons/proton_pure_eta"), eta);
        registryData.fill(HIST("data/pure/protons/proton_pure_dcaxy"), dcaxy);
        registryData.fill(HIST("data/pure/protons/proton_pure_dcaz"), dcaz);

        if (charge > 0) {
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_tpc"), pt, track.tpcNSigmaPr());
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_tof"), pt, track.tofNSigmaPr());
          }
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_pt"), pt);
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_eta"), eta);
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_dcaxy"), dcaxy);
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_dcaz"), dcaz);
        } else {
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_tpc"), pt, track.tpcNSigmaPr());
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_tof"), pt, track.tofNSigmaPr());
          }
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_pt"), pt);
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_eta"), eta);
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_dcaxy"), dcaxy);
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_dcaz"), dcaz);
        }
      }
    }
  }
  PROCESS_SWITCH(JetHadronsPid, processPureTracks, "Pure Tracks Analysis", true);

  void processJets(JetEvents::iterator const& collision,
                   soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                   soa::Join<aod::JetTracks, aod::JTrackPIs> const&,
                   HadronTracks const& globalTracks)
  {
    registryData.fill(HIST("data/jets/n_events_raw"), 1);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    float zVertex = collision.posZ();
    if (std::abs(zVertex) > zVtx) {
      return;
    }

    registryData.fill(HIST("data/jets/n_events"), 1);
    registryData.fill(HIST("data/jets/z_vtx"), zVertex);

    double centralRho = collision.rho();
    int jetsInCollision = 0;

    auto collTracks = globalTracks.sliceBy(tracksPerCollision, collision.collisionId());
    int tracksInJetsCollision = 0;

    for (auto const& jet : jets) {

      if (!isppRefAnalysis && ((std::abs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))) {
        continue;
      }
      if (isppRefAnalysis && std::abs(jet.eta()) > cfgEtaJetMax) {
        continue;
      }

      if (isppRefAnalysis && (jet.pt() < minJetPt || jet.pt() > maxJetPt)) {
        continue;
      }

      double ptSub = jet.pt() - (centralRho * jet.area());
      if (ptSub < 0) {
        ptSub = 0.0;
      }

      if (!isppRefAnalysis && (ptSub < minJetPt || ptSub > maxJetPt)) {
        continue;
      }

      double normalizedJetArea = jet.area() / (PI * rJet * rJet);
      if (applyAreaCut && normalizedJetArea < minNormalizedJetArea) {
        continue;
      }

      jetsInCollision++;

      registryData.fill(HIST("data/jets/jet_pt_subtracted"), ptSub);
      registryData.fill(HIST("data/jets/jet_pt_raw_vs_sub"), jet.pt(), ptSub);

      registryData.fill(HIST("data/jets/jet_eta"), jet.eta());
      registryData.fill(HIST("data/jets/jet_phi"), jet.phi());
      registryData.fill(HIST("data/jets/jet_area"), jet.area());
      registryData.fill(HIST("data/jets/jet_pt"), jet.pt());

      const double magnitudeThreshold = 1e-9;
      TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
      TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
      getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);

      if (ueAxis1.Mag() < magnitudeThreshold || ueAxis2.Mag() < magnitudeThreshold) {
        continue;
      }

      int selectedConstituentCount = 0;

      int nPionsInJet = 0;
      int nKaonsInJet = 0;
      int nProtonsInJet = 0;

      int nPionsInJetPos = 0;
      int nKaonsInJetPos = 0;
      int nProtonsInJetPos = 0;

      int nPionsInJetNeg = 0;
      int nKaonsInJetNeg = 0;
      int nProtonsInJetNeg = 0;

      for (auto const& jtrack : jet.tracks_as<soa::Join<aod::JetTracks, aod::JTrackPIs>>()) {

        auto track = jtrack.track_as<HadronTracks>();

        if (!passedTrackSelection(track)) {
          continue;
        }
        tracksInJetsCollision++;

        selectedConstituentCount++;
        const double p = track.p();

        registryData.fill(HIST("data/jets/tpc_signal_vs_p"), p, track.tpcSignal());

        if (track.hasTOF()) {
          registryData.fill(HIST("data/jets/tof_beta_vs_p"), p, track.beta());
        }

        double pt = track.pt();
        double eta = track.eta();
        double dcaxy = track.dcaXY();
        double dcaz = track.dcaZ();
        int charge = track.sign();
        double phi = track.phi();

        double deltaEtaJet = track.eta() - jet.eta();
        double deltaPhiJet = RecoDecay::constrainAngle(track.phi() - jet.phi(), -PI);

        registryData.fill(HIST("data/jets/cone/all_tracks"), deltaEtaJet, deltaPhiJet, pt);
        registryData.fill(HIST("data/jets/cone/all_tracks_2d"), deltaEtaJet, deltaPhiJet);

        PidResult pid = getPid(track);

        if (pid.isPion) {
          registryData.fill(HIST("data/jets/pions/pion_tpc_signal_vs_p"), track.p(), track.tpcSignal());

          if (track.hasTOF()) {registryData.fill(HIST("data/jets/pions/pion_tof_beta_vs_p"), p, track.beta());
          }

          if (charge > 0) {
            nPionsInJetPos++;

            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_tpc"), pt, track.tpcNSigmaPi());
            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/cone/pions/pion_pos"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/pions/pion_pos_2d"), deltaEtaJet, deltaPhiJet);

            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_tof"), pt, track.tofNSigmaPi());
            }
            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_pt"), pt);
            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_eta"), eta);
            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_dcaxy"), dcaxy);
            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_dcaz"), dcaz);
          } else {
            nPionsInJetNeg++;

            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_tpc"), pt, track.tpcNSigmaPi());
            registryData.fill(HIST("data/jets/cone/pions/pion_neg"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/pions/pion_neg_2d"), deltaEtaJet, deltaPhiJet);
            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_tof"), pt, track.tofNSigmaPi());
            }
            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_pt"), pt);
            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_eta"), eta);
            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_dcaxy"), dcaxy);
            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_dcaz"), dcaz);
          }
        }

        if (pid.isKaon) {
          registryData.fill(HIST("data/jets/kaons/kaon_tpc_signal_vs_p"), track.p(), track.tpcSignal());

          if (track.hasTOF()) {
            registryData.fill(HIST("data/jets/kaons/kaon_tof_beta_vs_p"), p, track.beta());
          }
          if (charge > 0) {
            nKaonsInJetPos++;

            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_tpc"), pt, track.tpcNSigmaKa());
            registryData.fill(HIST("data/jets/cone/kaons/kaon_pos"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/kaons/kaon_pos_2d"), deltaEtaJet, deltaPhiJet);

            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_tof"), pt, track.tofNSigmaKa());
            }
            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_pt"), pt);
            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_eta"), eta);
            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_dcaxy"), dcaxy);
            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_dcaz"), dcaz);
          } else {
            nKaonsInJetNeg++;

            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_tpc"), pt, track.tpcNSigmaKa());
            registryData.fill(HIST("data/jets/cone/kaons/kaon_neg"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/kaons/kaon_neg_2d"), deltaEtaJet, deltaPhiJet);
            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_tof"), pt, track.tofNSigmaKa());
            }
            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_pt"), pt);
            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_eta"), eta);
            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_dcaxy"), dcaxy);
            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_dcaz"), dcaz);
          }
        }

        if (pid.isProton) {
          registryData.fill(HIST("data/jets/protons/proton_tpc_signal_vs_p"), track.p(), track.tpcSignal());

          if (track.hasTOF()) {
            registryData.fill(
              HIST("data/jets/protons/proton_tof_beta_vs_p"),
              p,
              track.beta());
          }
          if (charge > 0) {
            nProtonsInJetPos++;

            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_tpc"), pt, track.tpcNSigmaPr());
            registryData.fill(HIST("data/jets/cone/protons/proton"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/protons/proton_2d"), deltaEtaJet, deltaPhiJet);
            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_tof"), pt, track.tofNSigmaPr());
            }
            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_pt"), pt);
            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_eta"), eta);
            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_dcaxy"), dcaxy);
            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_dcaz"), dcaz);
          } else {
            nProtonsInJetNeg++;

            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_tpc"), pt, track.tpcNSigmaPr());
            registryData.fill(HIST("data/jets/cone/protons/antiproton"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/protons/antiproton_2d"), deltaEtaJet, deltaPhiJet);
            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_tof"), pt, track.tofNSigmaPr());
            }
            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_pt"), pt);
            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_eta"), eta);
            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_dcaxy"), dcaxy);
            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_dcaz"), dcaz);
          }
        }
      }
      registryData.fill(HIST("data/jets/jet_n_constituents"), selectedConstituentCount);

      if (nPionsInJetPos > 0) {
        registryData.fill(HIST("data/jets/pions/pos/n_pions_per_jet"), nPionsInJetPos);
      }
      if (nKaonsInJetPos > 0) {
        registryData.fill(HIST("data/jets/kaons/pos/n_kaons_per_jet"), nKaonsInJetPos);
      }
      if (nProtonsInJetPos > 0) {
        registryData.fill(HIST("data/jets/protons/pos/n_protons_per_jet"), nProtonsInJetPos);
      }

      if (nPionsInJetNeg > 0) {
        registryData.fill(HIST("data/jets/pions/neg/n_pions_per_jet"), nPionsInJetNeg);
      }
      if (nKaonsInJetNeg > 0) {
        registryData.fill(HIST("data/jets/kaons/neg/n_kaons_per_jet"), nKaonsInJetNeg);
      }
      if (nProtonsInJetNeg > 0) {
        registryData.fill(HIST("data/jets/protons/neg/n_protons_per_jet"), nProtonsInJetNeg);
      }

      nPionsInJet = nPionsInJetPos + nPionsInJetNeg;
      nProtonsInJet = nProtonsInJetNeg + nProtonsInJetPos;
      nKaonsInJet = nKaonsInJetNeg + nKaonsInJetPos;

      registryData.fill(HIST("data/jets/composition/particle_multiplicity_per_jet"), 1.0, nPionsInJet);
      registryData.fill(HIST("data/jets/composition/particle_multiplicity_per_jet"), 2.0, nKaonsInJet);
      registryData.fill(HIST("data/jets/composition/particle_multiplicity_per_jet"), 3.0, nProtonsInJet);

      registryData.fill(HIST("data/jets/composition/nKaons_vs_nProtons"), nKaonsInJet, nProtonsInJet);
      registryData.fill(HIST("data/jets/composition/nPions_vs_nKaons"), nPionsInJet, nKaonsInJet);
      registryData.fill(HIST("data/jets/composition/nPions_vs_nProtons"), nPionsInJet, nProtonsInJet);
      registryData.fill(HIST("data/jets/composition/nPions_nKaons_nProtons"), nPionsInJet, nKaonsInJet, nProtonsInJet);

      switch (nProtonsInJet) {
        case 1:
          registryData.fill(HIST("data/jets/composition/nKaons_vs_nPions_1proton"), nKaonsInJet, nPionsInJet);
          break;

        case 2:
          registryData.fill(HIST("data/jets/composition/nKaons_vs_nPions_2protons"), nKaonsInJet, nPionsInJet);
          break;

        case 3:
          registryData.fill(HIST("data/jets/composition/nKaons_vs_nPions_3protons"), nKaonsInJet, nPionsInJet);
          break;

        case 4:
          registryData.fill(HIST("data/jets/composition/nKaons_vs_nPions_4protons"), nKaonsInJet, nPionsInJet);
          break;

        default:
          break;
      }

      int nTracksOut = 0;

      for (auto const& track : collTracks) {

        if (!passedTrackSelection(track)) {
          continue;
        }

        double deltaEtaUe1 = track.eta() - ueAxis1.Eta();
        double deltaPhiUe1 = RecoDecay::constrainAngle(track.phi() - ueAxis1.Phi(), -PI);
        double deltaRUe1 = std::hypot(deltaEtaUe1, deltaPhiUe1);

        double deltaEtaUe2 = track.eta() - ueAxis2.Eta();
        double deltaPhiUe2 = RecoDecay::constrainAngle(track.phi() - ueAxis2.Phi(), -PI);
        double deltaRUe2 = std::hypot(deltaEtaUe2, deltaPhiUe2);

        if (deltaRUe1 > rJet && deltaRUe2 > rJet) {
          continue;
        }

        double pt = track.pt();
        double eta = track.eta();
        double dcaxy = track.dcaXY();
        double dcaz = track.dcaZ();
        int charge = track.sign();

        PidResult pid = getPid(track);
        nTracksOut++;

        if (pid.isPion) {
          registryData.fill(HIST("data/ue/pions/pion_ue_tpc"), pt, track.tpcNSigmaPi());
          if (track.hasTOF()) {
            registryData.fill(HIST("data/ue/pions/pion_ue_tof"), pt, track.tofNSigmaPi());
          }

          registryData.fill(HIST("data/ue/pions/pion_ue_pt"), pt);
          registryData.fill(HIST("data/ue/pions/pion_ue_eta"), eta);
          registryData.fill(HIST("data/ue/pions/pion_ue_dcaxy"), dcaxy);
          registryData.fill(HIST("data/ue/pions/pion_ue_dcaz"), dcaz);

          if (charge > 0) {
            registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_tpc"), pt, track.tpcNSigmaPi());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_tof"), pt, track.tofNSigmaPi());
            }
            registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_pt"), pt);
            registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_eta"), eta);
            registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_dcaxy"), dcaxy);
            registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_dcaz"), dcaz);
          } else {
            registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_tpc"), pt, track.tpcNSigmaPi());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_tof"), pt, track.tofNSigmaPi());
            }
            registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_pt"), pt);
            registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_eta"), eta);
            registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_dcaxy"), dcaxy);
            registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_dcaz"), dcaz);
          }
        }

        if (pid.isKaon) {
          registryData.fill(HIST("data/ue/kaons/kaon_ue_tpc"), pt, track.tpcNSigmaKa());
          if (track.hasTOF()) {
            registryData.fill(HIST("data/ue/kaons/kaon_ue_tof"), pt, track.tofNSigmaKa());
          }

          registryData.fill(HIST("data/ue/kaons/kaon_ue_pt"), pt);
          registryData.fill(HIST("data/ue/kaons/kaon_ue_eta"), eta);
          registryData.fill(HIST("data/ue/kaons/kaon_ue_dcaxy"), dcaxy);
          registryData.fill(HIST("data/ue/kaons/kaon_ue_dcaz"), dcaz);

          if (charge > 0) {
            registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_tpc"), pt, track.tpcNSigmaKa());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_tof"), pt, track.tofNSigmaKa());
            }
            registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_pt"), pt);
            registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_eta"), eta);
            registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_dcaxy"), dcaxy);
            registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_dcaz"), dcaz);
          } else {
            registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_tpc"), pt, track.tpcNSigmaKa());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_tof"), pt, track.tofNSigmaKa());
            }
            registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_pt"), pt);
            registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_eta"), eta);
            registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_dcaxy"), dcaxy);
            registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_dcaz"), dcaz);
          }
        }

        if (pid.isProton) {
          registryData.fill(HIST("data/ue/protons/proton_ue_tpc"), pt, track.tpcNSigmaPr());
          if (track.hasTOF()) {
            registryData.fill(HIST("data/ue/protons/proton_ue_tof"), pt, track.tofNSigmaPr());
          }

          registryData.fill(HIST("data/ue/protons/proton_ue_pt"), pt);
          registryData.fill(HIST("data/ue/protons/proton_ue_eta"), eta);
          registryData.fill(HIST("data/ue/protons/proton_ue_dcaxy"), dcaxy);
          registryData.fill(HIST("data/ue/protons/proton_ue_dcaz"), dcaz);

          if (charge > 0) {
            registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_tpc"), pt, track.tpcNSigmaPr());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_tof"), pt, track.tofNSigmaPr());
            }
            registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_pt"), pt);
            registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_eta"), eta);
            registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_dcaxy"), dcaxy);
            registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_dcaz"), dcaz);
          } else {
            registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_tpc"), pt, track.tpcNSigmaPr());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_tof"), pt, track.tofNSigmaPr());
            }
            registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_pt"), pt);
            registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_eta"), eta);
            registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_dcaxy"), dcaxy);
            registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_dcaz"), dcaz);
          }
        }
      }
      registryData.fill(HIST("data/ue/tracks_n_in_ue"), nTracksOut);
    }
    registryData.fill(HIST("data/jets/collision_multiplicity"), jetsInCollision);
  }
  PROCESS_SWITCH(JetHadronsPid, processJets, "Jets Analysis", true);

  void processJetsMCD(MCDJetEvents::iterator const& collision, aod::JetMcCollisions const&, ChargedMCDJets const& detectorLevelJets, soa::Join<aod::JetTracks, aod::JTrackPIs> const&, HadronTracksMC const& globalTracks, aod::McParticles const& mcParticles)
{
  if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
    return;
  }

  if (std::abs(collision.posZ()) > zVtx) {
    return;
  }

  auto collTracks = globalTracks.sliceBy(mcTracksPerCollision, collision.collisionId());
  auto collMcParticles = mcParticles.sliceBy(mcParticlesPerCollision, collision.collisionId());

  registryData.fill(HIST("mc/jets/detector/n_events"), 1);

  for (auto const& jet : detectorLevelJets) {

    if (!isppRefAnalysis && ((std::abs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))) {
      continue;
    }

    if (isppRefAnalysis && std::abs(jet.eta()) > cfgEtaJetMax) {
      continue;
    }

    if (jet.pt() < minJetPt || jet.pt() > maxJetPt) {
      continue;
    }

    const double normalizedJetArea = jet.area() / (PI * rJet * rJet);

    if (applyAreaCut && normalizedJetArea < minNormalizedJetArea) {
      continue;
    }

    registryData.fill(HIST("mc/jets/detector/jet_pt"), jet.pt());
    registryData.fill(HIST("mc/jets/detector/jet_eta"), jet.eta());
    registryData.fill(HIST("mc/jets/detector/jet_phi"), jet.phi());
    registryData.fill(HIST("mc/jets/detector/jet_area"), jet.area());
    registryData.fill(HIST("mc/jets/detector/jet_n_constituents"), static_cast<double>(jet.tracksIds().size()));

    TVector3 jetAxis(jet.px(), jet.py(), jet.pz());

    TVector3 ueAxis1(0, 0, 0);
    TVector3 ueAxis2(0, 0, 0);

    getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);

    int nSelectedConstituents = 0;

    // 0 = unknown, 1 = quark, 2 = gluon
    int jetFlavor = 0;
    double maxPartonPt = -1.0;

    for (auto const& part : collMcParticles) {

      const int pdg = std::abs(part.pdgCode());

      const bool isQuark = pdg >= 1 && pdg <= 6;

      const bool isGluon = pdg == 21;

      if (!isQuark && !isGluon) {
        continue;
      }

      const double dEta = part.eta() - jet.eta();

      const double dPhi = RecoDecay::constrainAngle(part.phi() - jet.phi(), -PI);

      const double dR = std::hypot(dEta, dPhi);

      if (dR < rJet && part.pt() > maxPartonPt) {
        maxPartonPt = part.pt();
        jetFlavor = isQuark ? 1 : 2;
      }
    }

    int nPionsInJetPos = 0;
    int nPionsInJetNeg = 0;

    int nKaonsInJetPos = 0;
    int nKaonsInJetNeg = 0;

    int nProtonsInJetPos = 0;
    int nProtonsInJetNeg = 0;

    for (auto const& jtrack : jet.tracks_as<soa::Join<aod::JetTracks, aod::JTrackPIs>>()) {

      auto track = jtrack.track_as<HadronTracksMC>();

      if (!passedTrackSelection(track)) {
        continue;
      }

      ++nSelectedConstituents;

      const double deltaEtaJet = track.eta() - jet.eta();

      const double deltaPhiJet = RecoDecay::constrainAngle(track.phi() - jet.phi(), -PI);

      const double pt = track.pt();
      const int charge = track.sign();

      registryData.fill(HIST("mc/jets/detector/cone/all_tracks_2d"), deltaEtaJet, deltaPhiJet);
      registryData.fill(HIST("mc/jets/detector/cone/all_tracks_3d"), deltaEtaJet, deltaPhiJet, pt);

      const PidResult pid = getPid(track);

      if (pid.isPion) {
        registryData.fill(HIST("mc/jets/detector/cone/pions_2d"), deltaEtaJet, deltaPhiJet);
        registryData.fill(HIST("mc/jets/detector/cone/pions_3d"), deltaEtaJet, deltaPhiJet, pt);
        registryData.fill(HIST("mc/jets/detector/pions/pion_pt"), pt);

        if (charge > 0) {
          ++nPionsInJetPos;

          registryData.fill(HIST("mc/jets/detector/pions/pos/pion_jet_pos_tpc"), pt, track.tpcNSigmaPi());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/jets/detector/pions/pos/pion_jet_pos_tof"), pt, track.tofNSigmaPi());
          }
          registryData.fill(HIST("mc/jets/detector/pions/pos/pion_jet_pos_pt"), pt);
          registryData.fill(HIST("mc/jets/detector/pions/pos/pion_jet_pos_eta"), track.eta());
          registryData.fill(HIST("mc/jets/detector/pions/pos/pion_jet_pos_dcaxy"), track.dcaXY());
          registryData.fill(HIST("mc/jets/detector/pions/pos/pion_jet_pos_dcaz"), track.dcaZ());
          registryData.fill(HIST("mc/jets/detector/pions/pos/pion_jet_pos_phi_vs_pt"), pt, track.phi());

        } else if (charge < 0) {
          ++nPionsInJetNeg;

          registryData.fill(HIST("mc/jets/detector/pions/neg/pion_jet_neg_tpc"), pt, track.tpcNSigmaPi());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/jets/detector/pions/neg/pion_jet_neg_tof"), pt, track.tofNSigmaPi());
          }
          registryData.fill(HIST("mc/jets/detector/pions/neg/pion_jet_neg_pt"), pt);
          registryData.fill(HIST("mc/jets/detector/pions/neg/pion_jet_neg_eta"), track.eta());
          registryData.fill(HIST("mc/jets/detector/pions/neg/pion_jet_neg_dcaxy"), track.dcaXY());
          registryData.fill(HIST("mc/jets/detector/pions/neg/pion_jet_neg_dcaz"), track.dcaZ());
          registryData.fill(HIST("mc/jets/detector/pions/neg/pion_jet_neg_phi_vs_pt"), pt, track.phi());
        }

      } else if (pid.isKaon) {
        registryData.fill(HIST("mc/jets/detector/cone/kaons_2d"), deltaEtaJet, deltaPhiJet);
        registryData.fill(HIST("mc/jets/detector/cone/kaons_3d"), deltaEtaJet, deltaPhiJet, pt);
        registryData.fill(HIST("mc/jets/detector/kaons/kaon_pt"), pt);

        if (charge > 0) {
          ++nKaonsInJetPos;

          registryData.fill(HIST("mc/jets/detector/kaons/pos/kaon_jet_pos_tpc"), pt, track.tpcNSigmaKa());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/jets/detector/kaons/pos/kaon_jet_pos_tof"), pt, track.tofNSigmaKa());
          }
          registryData.fill(HIST("mc/jets/detector/kaons/pos/kaon_jet_pos_pt"), pt);
          registryData.fill(HIST("mc/jets/detector/kaons/pos/kaon_jet_pos_eta"), track.eta());
          registryData.fill(HIST("mc/jets/detector/kaons/pos/kaon_jet_pos_dcaxy"), track.dcaXY());
          registryData.fill(HIST("mc/jets/detector/kaons/pos/kaon_jet_pos_dcaz"), track.dcaZ());
          registryData.fill(HIST("mc/jets/detector/kaons/pos/kaon_jet_pos_phi_vs_pt"), pt, track.phi());

        } else if (charge < 0) {
          ++nKaonsInJetNeg;

          registryData.fill(HIST("mc/jets/detector/kaons/neg/kaon_jet_neg_tpc"), pt, track.tpcNSigmaKa());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/jets/detector/kaons/neg/kaon_jet_neg_tof"), pt, track.tofNSigmaKa());
          }
          registryData.fill(HIST("mc/jets/detector/kaons/neg/kaon_jet_neg_pt"), pt);
          registryData.fill(HIST("mc/jets/detector/kaons/neg/kaon_jet_neg_eta"), track.eta());
          registryData.fill(HIST("mc/jets/detector/kaons/neg/kaon_jet_neg_dcaxy"), track.dcaXY());
          registryData.fill(HIST("mc/jets/detector/kaons/neg/kaon_jet_neg_dcaz"), track.dcaZ());
          registryData.fill(HIST("mc/jets/detector/kaons/neg/kaon_jet_neg_phi_vs_pt"), pt, track.phi());
        }

      } else if (pid.isProton) {
        registryData.fill(HIST("mc/jets/detector/cone/protons_2d"), deltaEtaJet, deltaPhiJet);
        registryData.fill(HIST("mc/jets/detector/cone/protons_3d"), deltaEtaJet, deltaPhiJet, pt);
        registryData.fill(HIST("mc/jets/detector/protons/proton_pt"), pt);

        if (charge > 0) {
          ++nProtonsInJetPos;

          registryData.fill(HIST("mc/jets/detector/protons/pos/proton_jet_pos_tpc"), pt, track.tpcNSigmaPr());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/jets/detector/protons/pos/proton_jet_pos_tof"), pt, track.tofNSigmaPr());
          }
          registryData.fill(HIST("mc/jets/detector/protons/pos/proton_jet_pos_pt"), pt);
          registryData.fill(HIST("mc/jets/detector/protons/pos/proton_jet_pos_eta"), track.eta());
          registryData.fill(HIST("mc/jets/detector/protons/pos/proton_jet_pos_dcaxy"), track.dcaXY());
          registryData.fill(HIST("mc/jets/detector/protons/pos/proton_jet_pos_dcaz"), track.dcaZ());
          registryData.fill(HIST("mc/jets/detector/protons/pos/proton_jet_pos_phi_vs_pt"), pt, track.phi());

        } else if (charge < 0) {
          ++nProtonsInJetNeg;

          registryData.fill(HIST("mc/jets/detector/protons/neg/antiproton_jet_neg_tpc"), pt, track.tpcNSigmaPr());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/jets/detector/protons/neg/antiproton_jet_neg_tof"), pt, track.tofNSigmaPr());
          }
          registryData.fill(HIST("mc/jets/detector/protons/neg/antiproton_jet_neg_pt"), pt);
          registryData.fill(HIST("mc/jets/detector/protons/neg/antiproton_jet_neg_eta"), track.eta());
          registryData.fill(HIST("mc/jets/detector/protons/neg/antiproton_jet_neg_dcaxy"), track.dcaXY());
          registryData.fill(HIST("mc/jets/detector/protons/neg/antiproton_jet_neg_dcaz"), track.dcaZ());
          registryData.fill(HIST("mc/jets/detector/protons/neg/antiproton_jet_neg_phi_vs_pt"), pt, track.phi());
        }
      }
    }

    const int nPionsInJet =
      nPionsInJetPos + nPionsInJetNeg;

    const int nKaonsInJet =
      nKaonsInJetPos + nKaonsInJetNeg;

    const int nProtonsInJet =
      nProtonsInJetPos + nProtonsInJetNeg;

    registryData.fill(HIST("mc/jets/detector/jet_n_selected_constituents"), nSelectedConstituents);
    registryData.fill(HIST("mc/jets/detector/composition/n_pions_per_jet"), nPionsInJet);
    registryData.fill(HIST("mc/jets/detector/composition/n_kaons_per_jet"), nKaonsInJet);
    registryData.fill(HIST("mc/jets/detector/composition/n_protons_per_jet"), nProtonsInJet);

    if (nPionsInJetPos > 0) {registryData.fill(HIST("mc/jets/detector/composition/pions/pos/n_pions_per_jet"), nPionsInJetPos);}
    if (nPionsInJetNeg > 0) {registryData.fill(HIST("mc/jets/detector/composition/pions/neg/n_pions_per_jet"), nPionsInJetNeg);}
    if (nKaonsInJetPos > 0) {registryData.fill(HIST("mc/jets/detector/composition/kaons/pos/n_kaons_per_jet"), nKaonsInJetPos);}
    if (nKaonsInJetNeg > 0) {registryData.fill(HIST("mc/jets/detector/composition/kaons/neg/n_kaons_per_jet"), nKaonsInJetNeg);}
    if (nProtonsInJetPos > 0 ) {registryData.fill(HIST("mc/jets/detector/composition/protons/pos/n_protons_per_jet"), nProtonsInJetPos);}
    if (nProtonsInJetNeg> 0) {registryData.fill(HIST("mc/jets/detector/composition/protons/neg/n_antiprotons_per_jet"), nProtonsInJetNeg);}

    if (jetFlavor == 1) {
      registryData.fill(HIST("mc/jets/detector/flavor/quark_jet_pt"), jet.pt());
      registryData.fill(HIST("mc/jets/detector/flavor/quark_pions_per_jet"), nPionsInJet);
      registryData.fill(HIST("mc/jets/detector/flavor/quark_kaons_per_jet"), nKaonsInJet);
      registryData.fill(HIST("mc/jets/detector/flavor/quark_protons_per_jet"), nProtonsInJet);

    } else if (jetFlavor == 2) {
      registryData.fill(HIST("mc/jets/detector/flavor/gluon_jet_pt"), jet.pt());
      registryData.fill(HIST("mc/jets/detector/flavor/gluon_pions_per_jet"), nPionsInJet);
      registryData.fill(HIST("mc/jets/detector/flavor/gluon_kaons_per_jet"), nKaonsInJet);
      registryData.fill(HIST("mc/jets/detector/flavor/gluon_protons_per_jet"), nProtonsInJet);
    }

    registryData.fill(HIST("mc/jets/detector/composition/n_pions_vs_n_kaons"), nPionsInJet, nKaonsInJet);
    registryData.fill(HIST("mc/jets/detector/composition/n_pions_vs_n_protons"), nPionsInJet, nProtonsInJet);
    registryData.fill(HIST("mc/jets/detector/composition/n_kaons_vs_n_protons"), nKaonsInJet, nProtonsInJet);
    registryData.fill(HIST("mc/jets/detector/composition/n_pions_n_kaons_n_protons"), nPionsInJet, nKaonsInJet, nProtonsInJet);

    int nTracksInUe = 0;
    int nPionsInUe = 0;
    int nKaonsInUe = 0;
    int nProtonsInUe = 0;

    for (auto const& track : collTracks) {

      if (!passedTrackSelection(track)) {
        continue;
      }

      const double deltaEtaUe1 = track.eta() - ueAxis1.Eta();
      const double deltaPhiUe1 = RecoDecay::constrainAngle(track.phi() - ueAxis1.Phi(), -PI);
      const double deltaRUe1 = std::hypot(deltaEtaUe1, deltaPhiUe1);
      const double deltaEtaUe2 = track.eta() - ueAxis2.Eta();
      const double deltaPhiUe2 = RecoDecay::constrainAngle(track.phi() - ueAxis2.Phi(), -PI);
      const double deltaRUe2 = std::hypot(deltaEtaUe2, deltaPhiUe2);
      const bool inUe1 = deltaRUe1 <= rJet;
      const bool inUe2 = deltaRUe2 <= rJet;

      if (!inUe1 && !inUe2) {
        continue;
      }

      const bool useUe1 = inUe1 && (!inUe2 || deltaRUe1 <= deltaRUe2);

      const double deltaEtaUe = useUe1 ? deltaEtaUe1 : deltaEtaUe2;
      const double deltaPhiUe = useUe1 ? deltaPhiUe1 : deltaPhiUe2;

      const double pt = track.pt();

      ++nTracksInUe;

      registryData.fill(HIST("mc/jets/detector/ue/all_tracks_2d"), deltaEtaUe, deltaPhiUe);

      const PidResult pid = getPid(track);

      if (pid.isPion) {

        ++nPionsInUe;

        registryData.fill(HIST("mc/jets/detector/ue/pions_2d"), deltaEtaUe, deltaPhiUe);
        registryData.fill(HIST("mc/jets/detector/ue/pion_pt"), pt);

      } else if (pid.isKaon) {

        ++nKaonsInUe;

        registryData.fill(HIST("mc/jets/detector/ue/kaons_2d"), deltaEtaUe, deltaPhiUe);
        registryData.fill(HIST("mc/jets/detector/ue/kaon_pt"), pt);

      } else if (pid.isProton) {

        ++nProtonsInUe;
        registryData.fill(HIST("mc/jets/detector/ue/protons_2d"), deltaEtaUe, deltaPhiUe);
        registryData.fill(HIST("mc/jets/detector/ue/proton_pt"), pt);
      }
    }
    registryData.fill(HIST("mc/jets/detector/ue/n_tracks_per_jet"), nTracksInUe);
    registryData.fill(HIST("mc/jets/detector/ue/n_pions_per_jet"), nPionsInUe);
    registryData.fill(HIST("mc/jets/detector/ue/n_kaons_per_jet"), nKaonsInUe);
    registryData.fill(HIST("mc/jets/detector/ue/n_protons_per_jet"), nProtonsInUe);
  }
}

  PROCESS_SWITCH(JetHadronsPid, processJetsMCD, "Detector-level charged jets in MC", false);


  void processMC(StandardEvents::iterator const& collision, HadronTracksMC const& tracks, aod::McParticles const&)
  {
    if (!collision.sel8() || std::abs(collision.posZ()) > zVtx) {
      return;
    }
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return;
    }

    registryData.fill(HIST("mc/reconstruction/collisions/n_events"), 1);
    registryData.fill(HIST("mc/reconstruction/collisions/z_vertex"), collision.posZ());
    registryData.fill(HIST("mc/reconstruction/collisions/x_vertex"), collision.posX());
    registryData.fill(HIST("mc/reconstruction/collisions/y_vertex"), collision.posY());
    registryData.fill(HIST("mc/reconstruction/collisions/xy_vertex"), collision.posX(), collision.posY());
    registryData.fill(HIST("mc/reconstruction/collisions/n_contributors"), collision.numContrib());

    int nSelectedTracks = 0;
    int nMcMatchedTracks = 0;

    for (auto const& track : tracks) {

      if (!passedTrackSelection(track)) {
        continue;
      }

      const double p = track.p();
      registryData.fill(HIST("mc/reconstruction/hadrons/tpc_signal_vs_p"), p, track.tpcSignal());

      if (track.hasTOF()) {
        registryData.fill(HIST("mc/reconstruction/hadrons/tof_beta_vs_p"), p, track.beta());
      }

      ++nSelectedTracks;

      double pt = track.pt();

      if (!track.has_mcParticle()) {
        continue;
      }

      ++nMcMatchedTracks;
      auto const& trueParticle = track.mcParticle();
      int pdg = std::abs(trueParticle.pdgCode());

      if (pdg == PDG_t::kPiPlus) {
        registryData.fill(HIST("mc/reconstruction/pions/tpc_signal_vs_p"), p,track.tpcSignal());

        if (track.hasTOF()) {
          registryData.fill(HIST("mc/reconstruction/pions/tof_beta_vs_p"), p, track.beta());
        }
      }

      if (pdg == PDG_t::kKPlus) {
        registryData.fill(
          HIST("mc/reconstruction/kaons/tpc_signal_vs_p"), p, track.tpcSignal());
        if (track.hasTOF()) {
          registryData.fill(HIST("mc/reconstruction/kaons/tof_beta_vs_p"), p, track.beta());
        }
      }

      if (pdg == PDG_t::kProton) {
        registryData.fill(HIST("mc/reconstruction/protons/tpc_signal_vs_p"), p, track.tpcSignal());

        if (track.hasTOF()) {
          registryData.fill(HIST("mc/reconstruction/protons/tof_beta_vs_p"), p, track.beta());
        }
      }

      int realpdg = trueParticle.pdgCode();
      bool isPrimary = trueParticle.isPhysicalPrimary();
      const int charge = track.sign();

      registryData.fill(HIST("mc/reconstruction/hadrons/rec_hadron_all"), pt);

      if (track.hasTOF()) {
        registryData.fill(HIST("mc/reconstruction/hadrons/rec_hadron_tof_match"), pt);
      }

      if (isPrimary) {
        registryData.fill(HIST("mc/reconstruction/hadrons/mc_rec_hadron_pt"), pt);
      } else {
        registryData.fill(HIST("mc/reconstruction/hadrons/mc_sec_hadron_pt"), pt);
      }

      PidResult pid = getPid(track);

      if (pid.isPion) {
        registryData.fill(HIST("mc/reconstruction/pions/rec_pion_all"), pt);
        registryData.fill(HIST("mc/reconstruction/pions/contamination_matrix_pion"), realpdg, pt);

        if (charge > 0) {

          registryData.fill(HIST("mc/reconstruction/pions/pos/rec_pion_pos_all"), pt);

          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/pions/pos/rec_pion_pos_tof_matched"), pt);
          }
          registryData.fill(HIST("mc/reconstruction/pions/pos/contamination_matrix_pion_pos"), realpdg, pt);

          registryData.fill(HIST("mc/reconstruction/pions/pos/rec_pion_pos_tpc"), pt, track.tpcNSigmaPi());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/pions/pos/rec_pion_pos_tof"), pt, track.tofNSigmaPi());
          }
        } else {
          registryData.fill(HIST("mc/reconstruction/pions/neg/rec_pion_neg_all"), pt);

          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/pions/neg/rec_pion_neg_tof_matched"), pt);
          }
          registryData.fill(HIST("mc/reconstruction/pions/neg/contamination_matrix_pion_neg"), realpdg, pt);

          registryData.fill(HIST("mc/reconstruction/pions/neg/rec_pion_neg_tpc"), pt, track.tpcNSigmaPi());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/pions/neg/rec_pion_neg_tof"), pt, track.tofNSigmaPi());
          }
        }

        if (isPrimary) {
          if (pdg == PDG_t::kPiPlus) {
            registryData.fill(HIST("mc/reconstruction/pions/mc_rec_pion_pt"), pt);
            if (charge > 0) {
              registryData.fill(HIST("mc/reconstruction/pions/pos/mc_rec_pion_pos_pt"), pt);
            } else {
              registryData.fill(HIST("mc/reconstruction/pions/neg/mc_rec_pion_neg_pt"), pt);
            }
          }
        } else {
          registryData.fill(HIST("mc/reconstruction/pions/mc_sec_pion_pt"), pt);
          if (charge > 0) {
            registryData.fill(HIST("mc/reconstruction/pions/pos/mc_sec_pion_pos_pt"), pt);
          } else {
            registryData.fill(HIST("mc/reconstruction/pions/neg/mc_sec_pion_neg_pt"), pt);
          }
        }
      }

      if (pid.isKaon) {
        registryData.fill(HIST("mc/reconstruction/kaons/rec_kaon_all"), pt);
        registryData.fill(HIST("mc/reconstruction/kaons/contamination_matrix_kaon"), realpdg, pt);

        if (charge > 0) {

          registryData.fill(HIST("mc/reconstruction/kaons/pos/rec_kaon_pos_all"), pt);

          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/kaons/pos/rec_kaon_pos_tof_matched"), pt);
          }
          registryData.fill(HIST("mc/reconstruction/kaons/pos/contamination_matrix_kaon_pos"), realpdg, pt);

          registryData.fill(HIST("mc/reconstruction/kaons/pos/rec_kaon_pos_tpc"), pt, track.tpcNSigmaKa());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/kaons/pos/rec_kaon_pos_tof"), pt, track.tofNSigmaKa());
          }
        } else {

          registryData.fill(HIST("mc/reconstruction/kaons/neg/rec_kaon_neg_all"), pt);

          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/kaons/neg/rec_kaon_neg_tof_matched"), pt);
          }
          registryData.fill(HIST("mc/reconstruction/kaons/neg/contamination_matrix_kaon_neg"), realpdg, pt);

          registryData.fill(HIST("mc/reconstruction/kaons/neg/rec_kaon_neg_tpc"), pt, track.tpcNSigmaKa());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/kaons/neg/rec_kaon_neg_tof"), pt, track.tofNSigmaKa());
          }
        }

        if (isPrimary) {
          if (pdg == PDG_t::kKPlus) {
            registryData.fill(HIST("mc/reconstruction/kaons/mc_rec_kaon_pt"), pt);
            if (charge > 0) {
              registryData.fill(HIST("mc/reconstruction/kaons/pos/mc_rec_kaon_pos_pt"), pt);
            } else {
              registryData.fill(HIST("mc/reconstruction/kaons/neg/mc_rec_kaon_neg_pt"), pt);
            }
          }
        } else {
          registryData.fill(HIST("mc/reconstruction/kaons/mc_sec_kaon_pt"), pt);
          if (charge > 0) {
            registryData.fill(HIST("mc/reconstruction/kaons/pos/mc_sec_kaon_pos_pt"), pt);
          } else {
            registryData.fill(HIST("mc/reconstruction/kaons/neg/mc_sec_kaon_neg_pt"), pt);
          }
        }
      }

      if (pid.isProton) {
        registryData.fill(HIST("mc/reconstruction/protons/rec_proton_all"), pt);
        registryData.fill(HIST("mc/reconstruction/protons/contamination_matrix_proton"), realpdg, pt);

        if (charge > 0) {
          registryData.fill(HIST("mc/reconstruction/protons/pos/rec_proton_pos_all"), pt);
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/protons/pos/rec_proton_pos_tof_matched"), pt);
          }
          registryData.fill(HIST("mc/reconstruction/protons/pos/contamination_matrix_proton_pos"), realpdg, pt);

          registryData.fill(HIST("mc/reconstruction/protons/pos/rec_proton_pos_tpc"), pt, track.tpcNSigmaPr());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/protons/pos/rec_proton_pos_tof"), pt, track.tofNSigmaPr());
          }
        } else {

          registryData.fill(HIST("mc/reconstruction/protons/neg/rec_proton_neg_all"), pt);
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/protons/neg/rec_proton_neg_tof_matched"), pt);
          }
          registryData.fill(HIST("mc/reconstruction/protons/neg/contamination_matrix_proton_neg"), realpdg, pt);

          registryData.fill(HIST("mc/reconstruction/protons/neg/rec_proton_neg_tpc"), pt, track.tpcNSigmaPr());
          if (track.hasTOF()) {
            registryData.fill(HIST("mc/reconstruction/protons/neg/rec_proton_neg_tof"), pt, track.tofNSigmaPr());
          }
        }

        if (isPrimary) {
          if (pdg == PDG_t::kProton) {
            registryData.fill(HIST("mc/reconstruction/protons/mc_rec_proton_pt"), pt);
            if (charge > 0) {
              registryData.fill(HIST("mc/reconstruction/protons/pos/mc_rec_proton_pos_pt"), pt);
            } else {
              registryData.fill(HIST("mc/reconstruction/protons/neg/mc_rec_proton_neg_pt"), pt);
            }
          }
        } else {
          registryData.fill(HIST("mc/reconstruction/protons/mc_sec_proton_pt"), pt);
          if (charge > 0) {
            registryData.fill(HIST("mc/reconstruction/protons/pos/mc_sec_proton_pos_pt"), pt);
          } else {
            registryData.fill(HIST("mc/reconstruction/protons/neg/mc_sec_proton_neg_pt"), pt);
          }
        }
      }
    }
    registryData.fill(HIST("mc/reconstruction/collisions/n_selected_tracks"), nSelectedTracks);
    registryData.fill(HIST("mc/reconstruction/collisions/n_mc_matched_tracks"), nMcMatchedTracks);
  }
  PROCESS_SWITCH(JetHadronsPid, processMC, "Run on Monte Carlo", false);

  void processMCTruth(aod::McCollisions::iterator const& mcCollision,
                    aod::McParticles const& mcParticles)
{
  if (std::abs(mcCollision.posZ()) > zVtx) {
    return;
  }

  for (auto const& mcpart : mcParticles) {

    if (!mcpart.isPhysicalPrimary()) {
      continue;
    }

    if (mcpart.eta() < minEta || mcpart.eta() > maxEta) {
      continue;
    }

    const int originalPdg = mcpart.pdgCode();
    const int pdg = std::abs(originalPdg);

    const double pt = mcpart.pt();
    const double eta = mcpart.eta();
    const double phi = RecoDecay::constrainAngle(mcpart.phi(), 0.0);

    bool selectedHadron = false;

    if (pdg == PDG_t::kPiPlus && pt >= cfg.minPtPion && pt <= cfg.maxPtPion) {

      selectedHadron = true;

      registryData.fill(HIST("mc/truth/pions/mc_gen_pion_pt"), pt);
      registryData.fill(HIST("mc/truth/pions/mc_gen_pion_eta"), eta);
      registryData.fill(HIST("mc/truth/pions/mc_gen_pion_phi"), phi);
      registryData.fill(HIST("mc/truth/pions/mc_gen_pion_eta_phi"), phi, eta);
      registryData.fill(HIST("mc/truth/pions/mc_gen_pion_phi_vs_pt"), pt, phi);

      if (originalPdg > 0) {

        registryData.fill(HIST("mc/truth/pions/pos/mc_gen_pion_pos_pt"), pt);
        registryData.fill(HIST("mc/truth/pions/pos/mc_gen_pion_pos_eta"), eta);
        registryData.fill(HIST("mc/truth/pions/pos/mc_gen_pion_pos_phi"), phi);
        registryData.fill(HIST("mc/truth/pions/pos/mc_gen_pion_pos_eta_phi"), phi, eta);
        registryData.fill(HIST("mc/truth/pions/pos/mc_gen_pion_pos_phi_vs_pt"), pt, phi);

      } else {

        registryData.fill(HIST("mc/truth/pions/neg/mc_gen_pion_neg_pt"), pt);
        registryData.fill(HIST("mc/truth/pions/neg/mc_gen_pion_neg_eta"), eta);
        registryData.fill(HIST("mc/truth/pions/neg/mc_gen_pion_neg_phi"), phi);
        registryData.fill(HIST("mc/truth/pions/neg/mc_gen_pion_neg_eta_phi"), phi, eta);
        registryData.fill(HIST("mc/truth/pions/neg/mc_gen_pion_neg_phi_vs_pt"), pt, phi);
      }

    } else if (pdg == PDG_t::kKPlus && pt >= cfg.minPtKaon && pt <= cfg.maxPtKaon) {

      selectedHadron = true;

      registryData.fill(HIST("mc/truth/kaons/mc_gen_kaon_pt"), pt);
      registryData.fill(HIST("mc/truth/kaons/mc_gen_kaon_eta"), eta);
      registryData.fill(HIST("mc/truth/kaons/mc_gen_kaon_phi"), phi);
      registryData.fill(HIST("mc/truth/kaons/mc_gen_kaon_eta_phi"), phi, eta);
      registryData.fill(HIST("mc/truth/kaons/mc_gen_kaon_phi_vs_pt"), pt, phi);

      if (originalPdg > 0) {

        registryData.fill(HIST("mc/truth/kaons/pos/mc_gen_kaon_pos_pt"), pt);
        registryData.fill(HIST("mc/truth/kaons/pos/mc_gen_kaon_pos_eta"), eta);
        registryData.fill(HIST("mc/truth/kaons/pos/mc_gen_kaon_pos_phi"), phi);
        registryData.fill(HIST("mc/truth/kaons/pos/mc_gen_kaon_pos_eta_phi"), phi, eta);
        registryData.fill(HIST("mc/truth/kaons/pos/mc_gen_kaon_pos_phi_vs_pt"), pt, phi);

      } else {

        registryData.fill(HIST("mc/truth/kaons/neg/mc_gen_kaon_neg_pt"), pt);
        registryData.fill(HIST("mc/truth/kaons/neg/mc_gen_kaon_neg_eta"), eta);
        registryData.fill(HIST("mc/truth/kaons/neg/mc_gen_kaon_neg_phi"), phi);
        registryData.fill(HIST("mc/truth/kaons/neg/mc_gen_kaon_neg_eta_phi"), phi, eta);
        registryData.fill(HIST("mc/truth/kaons/neg/mc_gen_kaon_neg_phi_vs_pt"), pt, phi);
      }

    } else if (pdg == PDG_t::kProton && pt >= cfg.minPtProton && pt <= cfg.maxPtProton) {

      selectedHadron = true;

      registryData.fill(HIST("mc/truth/protons/mc_gen_proton_pt"), pt);
      registryData.fill(HIST("mc/truth/protons/mc_gen_proton_eta"), eta);
      registryData.fill(HIST("mc/truth/protons/mc_gen_proton_phi"), phi);
      registryData.fill(HIST("mc/truth/protons/mc_gen_proton_eta_phi"), phi, eta);
      registryData.fill(HIST("mc/truth/protons/mc_gen_proton_phi_vs_pt"), pt, phi);

      if (originalPdg > 0) {

        registryData.fill(HIST("mc/truth/protons/pos/mc_gen_proton_pos_pt"), pt);
        registryData.fill(HIST("mc/truth/protons/pos/mc_gen_proton_pos_eta"), eta);
        registryData.fill(HIST("mc/truth/protons/pos/mc_gen_proton_pos_phi"), phi);
        registryData.fill(HIST("mc/truth/protons/pos/mc_gen_proton_pos_eta_phi"), phi, eta);
        registryData.fill(HIST("mc/truth/protons/pos/mc_gen_proton_pos_phi_vs_pt"), pt, phi);

      } else {

        registryData.fill(HIST("mc/truth/protons/neg/mc_gen_proton_neg_pt"), pt);
        registryData.fill(HIST("mc/truth/protons/neg/mc_gen_proton_neg_eta"), eta);
        registryData.fill(HIST("mc/truth/protons/neg/mc_gen_proton_neg_phi"), phi);
        registryData.fill(HIST("mc/truth/protons/neg/mc_gen_proton_neg_eta_phi"), phi, eta);
        registryData.fill(HIST("mc/truth/protons/neg/mc_gen_proton_neg_phi_vs_pt"), pt, phi);
      }
    }

    if (selectedHadron) {
      registryData.fill(HIST("mc/truth/hadrons/mc_gen_hadron_pt"), pt);
      registryData.fill(HIST("mc/truth/hadrons/mc_gen_hadron_eta"), eta);
      registryData.fill(HIST("mc/truth/hadrons/mc_gen_hadron_phi"), phi);
      registryData.fill(HIST("mc/truth/hadrons/mc_gen_hadron_eta_phi"), phi, eta);
      registryData.fill(HIST("mc/truth/hadrons/mc_gen_hadron_phi_vs_pt"), pt, phi);
    }
  }
}

  PROCESS_SWITCH(JetHadronsPid, processMCTruth, "Run on Monte Carlo (Pure Truth)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetHadronsPid>(cfgc)};
}
