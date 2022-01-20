function IDL = importIDL(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  IDL = importIDL("G:\Shared drives\Science\IDL\4. Measures (SOP, questionnaires, psi tasks, data)\2. Questionnaires\Survey Data\Main Data\Master IDL Database 5_27_20.xlsx", "master sheet", [2, 1837]);

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, 1837];
end

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 368);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":ND" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["IDL_ID", "IDL_Sess", "RespodentID", "CollectorID", "StartDate_PST", "EndDate_PST", "IPAddress", "EmailAddress", "FirstName", "LastName", "CustomData1", "Lang", "Age", "Gender", "Country", "state", "city", "Att_1", "Md_st", "Md_ar", "Workshop", "FLFEserv", "UseMo", "UseYr", "FLFEBeta", "SrtedWrkshp", "Brght_PG", "Brght_S", "Brght_B", "Brght_A", "Brght_RR", "Brght_txt", "Ed", "Rel", "Race", "Race_other", "Income", "Housenum", "Setting", "Slp_Mnth", "Pn_Mnth", "Dx_MD", "Dx_anx", "DX_BD", "Dx_mania", "DX_psychschiz", "Dx_add", "Dx_PTSD", "Dx_DNA", "DX_txt", "Psych", "GnHlth", "FllwdHealth", "HlthCondDes_txt", "PI1", "PI2", "PI3", "PI4", "PI5", "PI6", "PI7", "PI8", "PI9", "PI10", "Tab1", "Tab1a", "Tab2", "Tab2a", "SpChld", "SpChld_Imp", "SpCurrent", "SpCurrent_Imp", "Cnt_PGWyr", "Cnt_PWGlf", "Meditate", "Med_minsess", "Med_often", "Med_yrs", "Int_SN", "Int_SO", "Caff", "Slp_24", "Pn_24", "AIOS", "Happ", "Enjoy", "Smile", "Stress", "Worry", "Sad", "C1", "C2", "C3", "C4", "C5", "Cr", "PANAS1", "PANAS2", "PANAS3", "PANAS4", "PANAS5", "PANAS6", "PANAS7", "PANAS8", "PANAS9", "PANAS10", "PB1", "PE1", "PB2", "PE2", "PB3", "PE3", "PB4", "PE4", "PB5", "PE5", "PB6", "PE6", "PB7", "PE7", "PB8", "PE8", "PB9", "PE9", "PB10", "PE10", "Synch1", "Int1", "Int2", "Int3", "Int4", "Int5", "Int6", "Int7", "Int8", "Int9", "Int10", "Int11", "Int12", "Int13", "Int14", "Int15", "Srvey_Length", "PosChg", "ClarIns", "CIMnExp_describe", "BehChg", "CxnTeach_2", "MnExp_wldan", "MnExpO", "MnExpPP", "CxnPpt", "DfcltChall", "CxnNat", "MysSpir", "MnExp_Text", "PE1_2", "PE2_2", "PE3_2", "PE4_2", "PE5_2", "PE6_2", "PE7_2", "PE8_2", "PE9_2", "PE10_2", "Synch_2", "tm4wrkshp", "Ease4Wrkshp", "Fdbk_txt", "Earthrise", "BrnHrt", "ExtraFld", "ForresearchassistantuseonlyCopypastethefullfilenamehereoncetheM", "Muse_Focus", "Audio", "Distr", "Userid", "alternateuse", "alternateuseimage", "session", "bubble_attention_mean", "bubble_attention_std", "bubble_baseline_mean", "bubble_baseline_std", "alternate_use_answer", "alternate_use_image", "bem_td_diff_log", "bem_td_diff_inv", "jar_count", "jar_image", "remote_viewing_hits", "time_estimation", "WorkshopLink", "jar_true", "Int", "jar_diff", "tm_diff", "Bubble_diff", "E", "A", "C", "N", "O", "PA", "NA", "Comp", "PB", "PBexp", "PANAS_P", "PANAS_N", "Image1", "Use1", "Use_Fluenc1", "Use_Flex1", "UseElab1", "Use_Org1", "WorkshopGroup", "Effct_INT", "Effct_Creat", "Effct_Trans", "Effct_Psi", "Effct_WB", "Effct_Other", "Frmt_Lec", "Frmt_SG", "Frmt_Pairs", "Frmt_Disc", "Frmt_Mvmt", "Frmt_Nat", "Frm_Other", "NP_Art", "NP_Death", "NP_Dream", "NP_EmbPrac", "NP_HH", "NP_Intention", "NP_Intuition", "NP_Meditation", "NP_Nature", "NP_Psi", "NP_PosPsych", "NP_Sound", "NP_Spirit", "NP_ASC", "NP_Tech", "NP_Other", "IDL_ID_Sess", "EEG_ID_Sess", "EEGid", "EEGsession", "EEGdate", "EEGtime", "badChan_1", "badChan_2", "badChan_3", "badChan_4", "lengthsec", "TP9_delta", "TP9_theta", "TP9_alpha", "TP9_beta", "TP9_gamma", "AF7_delta", "AF7_theta", "AF7_alpha", "AF7_beta", "AF7_gamma", "AF8_delta", "AF8_theta", "AF8_alpha", "AF8_beta", "AF8_gamma", "TP10_delta", "TP10_theta", "TP10_alpha", "TP10_beta", "TP10_gamma", "ECGID", "ECGsession", "ECGFileName", "ECGDate", "PRMDetrending", "PRMInterpRate", "PRMMinMaxHR", "PRMNNxxThreshold", "PRMVLFband", "PRMLFband", "PRMHFband", "PRMFreqPoints", "PRMFFTorLomb", "PRMWelchWindow", "PRMLombWindow", "PRMARspectrum", "PRMEntropy", "PRMDFAshortterm", "PRMDFAlongterm", "PRMRecurrencePlot", "PRMNbrSamples", "PRMArtifactCorrection", "OnsetOffset", "Artifact", "PNSindex", "SNSindex", "Stressindex", "MeanRRms", "SDNNms", "MeanHRbpm", "SDHRbpm", "MinHRbpm", "MaxHRbpm", "RMSSDms", "NNxxbeats", "pNNxx", "HRVtriangularindex", "TINNms", "SDANNms", "SDNNIms", "VLFpeak_ARHz", "LFpeak_ARHz", "HFpeak_ARHz", "VLFpow_ARms2", "LFpow_ARms2", "HFpow_ARms2", "VLFpow_ARlog", "LFpow_ARlog", "HFpow_ARlog", "VLFpow_AR", "LFpow_AR", "HFpow_AR", "LFpow_ARnu", "HFpow_ARnu", "TOTpow_ARms2", "LF_HF_ratio_AR", "EDRHz", "SD1ms", "SD2ms", "SD2_SD1_ratio", "ApEn", "SampEn", "D2", "DFA1", "DFA2", "RP_Lmeanbeats", "RP_Lmaxbeats", "RP_REC", "RP_DET", "RP_ShanEn", "MSE_1", "MSE_2", "MSE_3", "MSE_4", "MSE_5", "MSE_6", "MSE_7", "MSE_8", "MSE_9", "MSE_10", "MSE_11", "MSE_12", "MSE_13", "MSE_14", "MSE_15", "MSE_16", "MSE_17", "MSE_18", "MSE_19", "MSE_20"];
opts.VariableTypes = ["categorical", "double", "char", "char", "datetime", "datetime", "char", "char", "char", "char", "char", "categorical", "double", "categorical", "categorical", "categorical", "char", "double", "char", "char", "categorical", "char", "char", "char", "char", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "char", "double", "categorical", "categorical", "char", "categorical", "double", "categorical", "double", "double", "categorical", "categorical", "char", "char", "char", "categorical", "categorical", "categorical", "char", "categorical", "categorical", "char", "char", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "categorical", "categorical", "categorical", "double", "double", "categorical", "char", "categorical", "char", "double", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "char", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "double", "categorical", "double", "char", "categorical", "categorical", "categorical", "categorical", "categorical", "char", "categorical", "categorical", "char", "categorical", "categorical", "char", "categorical", "categorical", "char", "categorical", "categorical", "char", "categorical", "char", "double", "double", "char", "categorical", "char", "char", "char", "double", "char", "char", "categorical", "char", "char", "double", "double", "double", "double", "double", "char", "categorical", "double", "double", "double", "categorical", "double", "double", "char", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "char", "char", "double", "char", "double", "double", "double", "double", "categorical", "double", "double", "double", "double", "double", "char", "double", "double", "double", "double", "double", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "char", "char", "categorical", "double", "categorical", "char", "categorical", "categorical", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char"];

% Specify variable properties
opts = setvaropts(opts, ["RespodentID", "CollectorID", "EmailAddress", "FirstName", "LastName", "CustomData1", "city", "Md_st", "Md_ar", "FLFEserv", "UseMo", "UseYr", "FLFEBeta", "Brght_txt", "Race_other", "DX_BD", "Dx_mania", "DX_psychschiz", "DX_txt", "FllwdHealth", "HlthCondDes_txt", "Med_minsess", "Med_yrs", "PANAS1", "PANAS2", "PANAS3", "PANAS4", "PANAS5", "PANAS6", "PANAS7", "PANAS8", "PANAS9", "PANAS10", "Synch1", "BehChg", "DfcltChall", "MnExp_Text", "PE3_2", "PE6_2", "PE9_2", "Synch_2", "Fdbk_txt", "BrnHrt", "ExtraFld", "ForresearchassistantuseonlyCopypastethefullfilenamehereoncetheM", "Audio", "Distr", "alternateuse", "alternateuseimage", "alternate_use_answer", "WorkshopLink", "PANAS_P", "PANAS_N", "Use1", "Effct_Other", "IDL_ID_Sess", "EEG_ID_Sess", "EEGtime", "ECGID", "ECGsession", "ECGFileName", "ECGDate", "PRMDetrending", "PRMInterpRate", "PRMMinMaxHR", "PRMNNxxThreshold", "PRMVLFband", "PRMLFband", "PRMHFband", "PRMFreqPoints", "PRMFFTorLomb", "PRMWelchWindow", "PRMLombWindow", "PRMARspectrum", "PRMEntropy", "PRMDFAshortterm", "PRMDFAlongterm", "PRMRecurrencePlot", "PRMNbrSamples", "PRMArtifactCorrection", "OnsetOffset", "Artifact", "PNSindex", "SNSindex", "Stressindex", "MeanRRms", "SDNNms", "MeanHRbpm", "SDHRbpm", "MinHRbpm", "MaxHRbpm", "RMSSDms", "NNxxbeats", "pNNxx", "HRVtriangularindex", "TINNms", "SDANNms", "SDNNIms", "VLFpeak_ARHz", "LFpeak_ARHz", "HFpeak_ARHz", "VLFpow_ARms2", "LFpow_ARms2", "HFpow_ARms2", "VLFpow_ARlog", "LFpow_ARlog", "HFpow_ARlog", "VLFpow_AR", "LFpow_AR", "HFpow_AR", "LFpow_ARnu", "HFpow_ARnu", "TOTpow_ARms2", "LF_HF_ratio_AR", "EDRHz", "SD1ms", "SD2ms", "SD2_SD1_ratio", "ApEn", "SampEn", "D2", "DFA1", "DFA2", "RP_Lmeanbeats", "RP_Lmaxbeats", "RP_REC", "RP_DET", "RP_ShanEn", "MSE_1", "MSE_2", "MSE_3", "MSE_4", "MSE_5", "MSE_6", "MSE_7", "MSE_8", "MSE_9", "MSE_10", "MSE_11", "MSE_12", "MSE_13", "MSE_14", "MSE_15", "MSE_16", "MSE_17", "MSE_18", "MSE_19", "MSE_20"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["IDL_ID", "RespodentID", "CollectorID", "EmailAddress", "FirstName", "LastName", "CustomData1", "Lang", "Gender", "Country", "state", "city", "Md_st", "Md_ar", "Workshop", "FLFEserv", "UseMo", "UseYr", "FLFEBeta", "SrtedWrkshp", "Brght_PG", "Brght_S", "Brght_B", "Brght_A", "Brght_RR", "Brght_txt", "Rel", "Race", "Race_other", "Income", "Setting", "Dx_MD", "Dx_anx", "DX_BD", "Dx_mania", "DX_psychschiz", "Dx_add", "Dx_PTSD", "Dx_DNA", "DX_txt", "Psych", "GnHlth", "FllwdHealth", "HlthCondDes_txt", "SpChld", "SpChld_Imp", "SpCurrent", "SpCurrent_Imp", "Meditate", "Med_minsess", "Med_often", "Med_yrs", "Caff", "PANAS1", "PANAS2", "PANAS3", "PANAS4", "PANAS5", "PANAS6", "PANAS7", "PANAS8", "PANAS9", "PANAS10", "Synch1", "Srvey_Length", "ClarIns", "BehChg", "CxnTeach_2", "MnExp_wldan", "MnExpO", "MnExpPP", "CxnPpt", "DfcltChall", "CxnNat", "MysSpir", "MnExp_Text", "PE1_2", "PE2_2", "PE3_2", "PE4_2", "PE5_2", "PE6_2", "PE7_2", "PE8_2", "PE9_2", "PE10_2", "Synch_2", "Fdbk_txt", "Earthrise", "BrnHrt", "ExtraFld", "ForresearchassistantuseonlyCopypastethefullfilenamehereoncetheM", "Audio", "Distr", "Userid", "alternateuse", "alternateuseimage", "alternate_use_answer", "alternate_use_image", "jar_image", "WorkshopLink", "PANAS_P", "PANAS_N", "Use1", "WorkshopGroup", "Effct_Other", "Frm_Other", "NP_Other", "IDL_ID_Sess", "EEG_ID_Sess", "EEGid", "EEGdate", "EEGtime", "badChan_1", "badChan_2", "badChan_3", "badChan_4", "ECGID", "ECGsession", "ECGFileName", "ECGDate", "PRMDetrending", "PRMInterpRate", "PRMMinMaxHR", "PRMNNxxThreshold", "PRMVLFband", "PRMLFband", "PRMHFband", "PRMFreqPoints", "PRMFFTorLomb", "PRMWelchWindow", "PRMLombWindow", "PRMARspectrum", "PRMEntropy", "PRMDFAshortterm", "PRMDFAlongterm", "PRMRecurrencePlot", "PRMNbrSamples", "PRMArtifactCorrection", "OnsetOffset", "Artifact", "PNSindex", "SNSindex", "Stressindex", "MeanRRms", "SDNNms", "MeanHRbpm", "SDHRbpm", "MinHRbpm", "MaxHRbpm", "RMSSDms", "NNxxbeats", "pNNxx", "HRVtriangularindex", "TINNms", "SDANNms", "SDNNIms", "VLFpeak_ARHz", "LFpeak_ARHz", "HFpeak_ARHz", "VLFpow_ARms2", "LFpow_ARms2", "HFpow_ARms2", "VLFpow_ARlog", "LFpow_ARlog", "HFpow_ARlog", "VLFpow_AR", "LFpow_AR", "HFpow_AR", "LFpow_ARnu", "HFpow_ARnu", "TOTpow_ARms2", "LF_HF_ratio_AR", "EDRHz", "SD1ms", "SD2ms", "SD2_SD1_ratio", "ApEn", "SampEn", "D2", "DFA1", "DFA2", "RP_Lmeanbeats", "RP_Lmaxbeats", "RP_REC", "RP_DET", "RP_ShanEn", "MSE_1", "MSE_2", "MSE_3", "MSE_4", "MSE_5", "MSE_6", "MSE_7", "MSE_8", "MSE_9", "MSE_10", "MSE_11", "MSE_12", "MSE_13", "MSE_14", "MSE_15", "MSE_16", "MSE_17", "MSE_18", "MSE_19", "MSE_20"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "StartDate_PST", "InputFormat", "");
opts = setvaropts(opts, "EndDate_PST", "InputFormat", "");

% Import the data
IDL = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":ND" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    IDL = [IDL; tb]; %#ok<AGROW>
end

end