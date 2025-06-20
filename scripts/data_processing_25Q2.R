library(tidyverse)
library(magrittr)
library(scam)
library(parallel)
source("scripts/utilities.R")

#----
# LOAD THE RAW DATA
#----
compound_annotations <- data.table::fread("data/25Q2/PRISMOncologyReferenceCompoundList.csv")
inst_meta <- data.table::fread("data/25Q2/PRISMOncologyReferenceInstMeta.csv")
analyte_meta <- data.table::fread("data/25Q2/PRISMOncologyReferenceAnalyteMeta.csv")
lmfi <- data.table::fread("data/25Q2/PRISMOncologyReferenceLMFI.csv")

# -----
# FILTER WELLS AND ANALYTES WITH LOW BEAD-COUNTS
# ESTIMATE THE NOISE FLOOR FOR EACH WELL & CORRECT THE LMFI VALUES
# -----

# Estimate the noise floor 
NoiseFloorEstimates <- lmfi %>% 
  dplyr::filter(COUNT >= 10, LMFI >= 4) %>% 
  dplyr::group_by(prism_replicate, pert_well) %>% 
  dplyr::arrange(LMFI) %>% 
  dplyr::mutate(rank = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(rank <= 10) %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::group_by(pert_well, prism_replicate) %>%
  dplyr::summarize(nf_nonexpected = mean(LMFI[is.na(pool_id)]), nf_expected = mean(LMFI[(rank <= 1) & (!is.na(pool_id))])) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(noise_floor = ifelse(is.finite(nf_nonexpected), nf_nonexpected, nf_expected)) 



# Remove low bead-count or too dim points and correct for the background flourescence.
lmfi <- lmfi %>%  
  dplyr::left_join(analyte_meta) %>% 
  dplyr::filter(!is.na(pool_id))%>% 
  dplyr::group_by(prism_replicate, pert_well) %>%
  dplyr::filter(quantile(COUNT, probs = 0.25, na.rm = T) >= 10) %>% 
  dplyr::ungroup() %>%
  dplyr::filter(COUNT >= 10, LMFI >= 4) %>% 
  dplyr::inner_join(NoiseFloorEstimates) %>% 
  dplyr::filter(LMFI > noise_floor) %>% 
  dplyr::mutate(LMFI.NF = log2(2^LMFI - 2^noise_floor)) %>% 
  dplyr::distinct(prism_replicate, pert_well, analyte_id, screen, cellset, COUNT, LMFI, LMFI.NF)



# -----
# FILTER NEGATIVE CONTROL WELLS THAT CONTROL BARCODES DOESN'T COVER AT LEAST 25% OF THE RANGE OF CELL LINE BARCODES OR
# AND ANY WELL THAT HAS SPEARMAN CORRELATION LESS THAN 0.5 WITH THE MEDIAN OF CONTROL BARCODES ACROSS NEGATIVE CONTROL WELLS.
# ----

CB.quantiles <-lmfi %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::group_by(prism_replicate, pert_well) %>%
  dplyr::arrange(LMFI.NF) %>% 
  dplyr::mutate(quant = (1:n()) / n()) %>% 
  dplyr::filter(pool_id == "CTLBC") %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::distinct(prism_replicate, pert_well, pert_type,
                  depmap_id, LMFI.NF, quant) %>% 
  dplyr::group_by(prism_replicate, pert_well) %>% 
  dplyr::mutate(range = max(quant, na.rm = T) - min(quant, na.rm = T),
                n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter((range >= 0.25) | (pert_type != "ctl_vehicle"), n > 5) 


CB.quantiles <- CB.quantiles %>% 
  dplyr::filter(pert_type == "ctl_vehicle") %>% 
  dplyr::group_by(prism_replicate, pert_well) %>% 
  dplyr::mutate(m = median(LMFI.NF, na.rm = T)) %>%
  dplyr::group_by(prism_replicate) %>% 
  dplyr::mutate(mm = median(m, na.rm = T)) %>% 
  dplyr::group_by(prism_replicate, depmap_id) %>% 
  dplyr::summarise(m.LMFI = median(LMFI.NF - m + mm, na.rm = T)) %>% 
  dplyr::ungroup() %>%
  dplyr::left_join(CB.quantiles) 


lmfi <- CB.quantiles %>% 
  dplyr::group_by(prism_replicate, pert_well, pert_type) %>%
  dplyr::summarise(concordance = cor(m.LMFI, LMFI.NF, use = "p", method = "spearman")) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(concordance >= 0.5) %>%
  dplyr::distinct(prism_replicate, pert_well) %>%
  dplyr::left_join(lmfi) 

rm(CB.quantiles)


# ------
# DROP PLATES THAT DOESN'T HAVE AT LEAST 16 NEGCON WELLS
# ------

lmfi <- lmfi %>%
  dplyr::distinct(prism_replicate, pert_well) %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::count(prism_replicate, pert_type) %>% 
  tidyr::pivot_wider(names_from = pert_type, values_from = n) %>% 
  dplyr::filter(trt_poscon >= 16, ctl_vehicle >= 16) %>% 
  dplyr::distinct(prism_replicate) %>% 
  dplyr::left_join(lmfi)


# ----
# COMPUTE THE REFERENCE VALUES FOR CONTROL BARCODES
# ----

REF = lmfi %>%  
  dplyr::semi_join(dplyr::filter(inst_meta, pert_type == "ctl_vehicle")) %>% 
  dplyr::semi_join(dplyr::filter(analyte_meta, pool_id == "CTLBC")) %>% 
  dplyr::group_by(prism_replicate, pert_well) %>% 
  dplyr::mutate(m = median(LMFI.NF, na.rm = T)) %>% 
  dplyr::group_by(prism_replicate) %>%  
  dplyr::mutate(m = m - median(m, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(LMFI.NF = LMFI.NF - m) %>% 
  dplyr::group_by(prism_replicate, analyte_id) %>%
  dplyr::summarise(LMFI.ref = median(LMFI.NF, na.rm = T)) %>% 
  dplyr::ungroup()

# ----
# NORMALIZE USING SPLINE FITS
# ----

lmfi %<>%
  dplyr::left_join(REF) %>%
  dplyr::mutate(profile_id = paste0(prism_replicate, "::", pert_well))

profiles = lmfi$profile_id %>% unique

lmfi.normalized = list()
ix = 1
for (profile in profiles) {
  print(ix)
  temp = lmfi %>%
    dplyr::filter(profile_id == profile)
  
  
  g <- tryCatch(
    scam(
      LMFI.ref ~ s(LMFI.NF, bs = "micx", k = 4),
      data = temp %>% dplyr::filter(is.finite(LMFI.ref)),
      optimizer = "nlm.fd"
    ),
    error = function(e)
      NULL
  )
  
  if (!is.null(g)) {
    p = predict(g, temp, se.fit = T)
    temp$LMFI.normalized = c(p$fit)
    temp$LMFI.normalized.se = c(p$se.fit)
  }
  
  lmfi.normalized[[ix]] = temp
  ix = ix + 1
}

lmfi.normalized %<>% dplyr::bind_rows()

rm(REF, temp, p, g, ix, profile, profiles)

# ----
# FILTERING OUTLIER NEGATIVE CONTROL POOLS AND ALIGN THE REPLICATES
# ----

outlier.nc.pools <- lmfi.normalized %>% 
  dplyr::left_join(inst_meta)%>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::filter(pert_type == "ctl_vehicle", pool_id != "CTLBC") %>%
  dplyr::group_by(prism_replicate, analyte_id) %>% 
  dplyr::mutate(m.LMFI.normalized = median(LMFI.normalized, na.rm = T)) %>% 
  dplyr::group_by(CompoundPlate, prism_replicate, pert_well, pool_id, cellset) %>% 
  dplyr::summarize(r = cor(LMFI.normalized, m.LMFI.normalized, use = "p", method = "s"),
                   mLFC = median(LMFI.normalized - m.LMFI.normalized, na.rm =T),
                   n = n()) %>%
  dplyr::group_by(prism_replicate, pool_id) %>% 
  dplyr::mutate(r.z = c(scale(r))) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(pass = (abs(mLFC) < 1) | (r.z > -2)) %>% 
  dplyr::group_by(prism_replicate, pert_well) %>% 
  dplyr::mutate(well.pass.ratio = mean(pass, na.rm = T)) %>% 
  dplyr::group_by(prism_replicate, pool_id) %>% 
  dplyr::mutate(pool.pass.ratio = mean(pass, na.rm = T)) %>% 
  dplyr::filter((well.pass.ratio < .75) | (pool.pass.ratio < .75) | (!pass)) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(prism_replicate, pert_well, pool_id)

ctrl.pool.corrections <- lmfi.normalized %>% 
  dplyr::left_join(inst_meta)%>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::filter(pert_type == "ctl_vehicle") %>% 
  dplyr::anti_join(outlier.nc.pools) %>% 
  dplyr::group_by(prism_replicate, pert_well, pool_id, cellset, pert_type) %>% 
  dplyr::summarise(mLMFI = median(LMFI.normalized, na.rm = T)) %>% 
  dplyr::group_by(prism_replicate, pert_well) %>% 
  dplyr::mutate(dLMFI = mLMFI - median(mLMFI[pool_id == "CTLBC"], na.rm = T)) %>% 
  dplyr::group_by(prism_replicate, pool_id, pert_type) %>% 
  dplyr::mutate(rLMFI = median(dLMFI, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(delta = rLMFI - dLMFI) 

lmfi.normalized <- lmfi.normalized %>% 
  dplyr::left_join(analyte_meta) %>%  
  dplyr::left_join(ctrl.pool.corrections %>%
                     dplyr::distinct(pert_well, prism_replicate, pool_id, delta)) %>% 
  dplyr::mutate(delta = ifelse(is.finite(delta), delta,0)) %>% 
  dplyr::mutate(LMFI.normalized.corrected =  LMFI.normalized + delta) %>% 
  dplyr::distinct(prism_replicate, pert_well, analyte_id, cellset, screen, COUNT, 
                  LMFI, LMFI.NF, LMFI.normalized, LMFI.normalized.corrected)  


# ----
# QC TABLE
# ----

QC = lmfi.normalized %>%
  dplyr::left_join(inst_meta) %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::filter(pert_type %in% c("ctl_vehicle", "trt_poscon"), pool_id != "CTLBC") %>%
  dplyr::group_by(screen, CompoundPlate, prism_replicate, cellset, analyte_id, pool_id, depmap_id) %>%
  dplyr::filter(is.finite(LMFI.normalized.corrected)) %>%  
  dplyr::summarise(
    error_rate =  min(PRROC::roc.curve(scores.class0 = LMFI.normalized.corrected,
                                       weights.class0 = pert_type == "ctl_vehicle",
                                       curve = TRUE)$curve[,1] + 1 -
                        PRROC::roc.curve(scores.class0 = LMFI.normalized.corrected,
                                         weights.class0 = pert_type == "ctl_vehicle",
                                         curve = TRUE )$curve[,2])/2,
    NC.median = median(LMFI.normalized.corrected[pert_type == "ctl_vehicle"], na.rm = T),
    NC.mad = mad(LMFI.normalized.corrected[pert_type == "ctl_vehicle"], na.rm = T),
    PC.median = median(LMFI.normalized.corrected[pert_type == "trt_poscon"], na.rm = T),
    PC.mad = mad( LMFI.normalized.corrected[pert_type == "trt_poscon"], na.rm = T),
    NC.n = sum(pert_type == "ctl_vehicle", na.rm = T),
    PC.n = sum(pert_type == "trt_poscon", na.rm = T)) %>%
  dplyr::mutate(DR = NC.median - PC.median,
                SSMD = DR / sqrt(NC.mad ^ 2 + PC.mad ^ 2)) %>%
  dplyr::mutate(PASS = (error_rate <= 0.05) & (DR > 2) & (NC.n >= 16) & (PC.n >= 16)) %>%
  dplyr::distinct() %>%
  dplyr::ungroup()

QC %>% 
  write_csv("data/25Q2/PRISMOncologyReferenceQCTable.csv")

# ----
# COMPUTE LOG-FOLD-CHANGES
# ----

LFC <- lmfi.normalized %>%  
  dplyr::inner_join(QC %>% dplyr::select(analyte_id, pool_id, prism_replicate, NC.median, PASS)) %>%
  dplyr::mutate(LFC = LMFI.normalized.corrected - NC.median) %>%
  dplyr::distinct(prism_replicate, pert_well, analyte_id, cellset, pool_id, screen, LFC, PASS) 

# ----
# CORRECT FOR POOL EFFECTS FOR TREATMENTS 
# ---- 

trt.pool.corrections <- LFC %>%
  dplyr::left_join(inst_meta) %>% 
  dplyr::filter(pert_type == "trt_cp") %>% 
  dplyr::group_by(prism_replicate, pert_well, pool_id, CompoundPlate, cellset, SampleID, pert_dose) %>% 
  dplyr::summarise(mLFC = median(LFC, na.rm = T)) %>% 
  dplyr::group_by(CompoundPlate, SampleID, pert_dose, pool_id, cellset) %>%
  dplyr::mutate(delta = median(mLFC, na.rm = T) - mLFC) %>%
  dplyr::ungroup()


LFC <- LFC %>% 
  dplyr::left_join(trt.pool.corrections %>%
                     dplyr::distinct(pert_well, prism_replicate, pool_id, delta)) %>% 
  dplyr::mutate(delta = ifelse(is.finite(delta), delta,0)) %>% 
  dplyr::mutate(LFC_corrected =  LFC + delta) %>% 
  dplyr::distinct(prism_replicate, pert_well, analyte_id, cellset, pool_id, screen, LFC, LFC_corrected, PASS)

LFC %>% 
  dplyr::distinct(prism_replicate, pert_well, analyte_id, cellset, screen, LFC, LFC_corrected, PASS) %>% 
  write_csv("data/25Q2/PRISMOncologyReferenceLFC.csv")

# ----
# COLLAPSE REPLICATES 
# ----

LFC.collapsed <- LFC %>% 
  dplyr::filter(PASS) %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::filter(pert_type == "trt_cp", !is.na(depmap_id), pool_id != "CTLBC", is.finite(LFC)) %>%  
  dplyr::group_by(CompoundPlate, SampleID, pert_dose, pert_dose_unit, cellset, pool_id, depmap_id) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::summarise(LFC = median(LFC_corrected, na.rm = T)) %>% 
  dplyr::ungroup()

# ----
# FLAG AND FILTER OUTLIER TREATMENT WELLS
# ----

dose_indices <- inst_meta %>% 
  dplyr::filter(pert_type == "trt_cp") %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose) %>% 
  dplyr::group_by(CompoundPlate, SampleID) %>% 
  dplyr::arrange(pert_dose) %>% 
  dplyr::mutate(dose_ix = 1:n()) %>%
  dplyr::ungroup()


monotonicity_flags <- dose_indices %>% 
  dplyr::left_join(dose_indices %>% 
                     dplyr::rename(pert_dose.prev = pert_dose) %>% 
                     dplyr::mutate(dose_ix = dose_ix + 1)) %>%   
  dplyr::left_join(dose_indices %>% 
                     dplyr::rename(pert_dose.next = pert_dose) %>% 
                     dplyr::mutate(dose_ix = dose_ix - 1)) %>%   
  dplyr::select(-dose_ix) %>% 
  dplyr::left_join(LFC.collapsed) %>%  
  dplyr::left_join(LFC.collapsed %>% 
                     dplyr::rename(pert_dose.prev = pert_dose,  LFC.prev = LFC)) %>%  
  dplyr::left_join(LFC.collapsed %>% 
                     dplyr::rename(pert_dose.next = pert_dose,  LFC.next = LFC)) %>%  
  dplyr::mutate(flag.down = LFC < -2, flag.up =  LFC > -1,
                flag.down.prev = ifelse(is.na(pert_dose.prev), FALSE, LFC.prev < -2), 
                flag.up.prev =  ifelse(is.na(pert_dose.prev), TRUE, LFC.prev > -1),
                flag.down.next = ifelse(is.na(pert_dose.next), TRUE, LFC.next < -2),
                flag.up.next =  ifelse(is.na(pert_dose.next), FALSE, LFC.next > -1)) %>% 
  dplyr::mutate(flag1 = flag.down & flag.up.next & flag.up.prev,
                flag2 = flag.up & flag.down.next & flag.down.prev)


outlier.trt.pools <- monotonicity_flags %>% 
  dplyr::group_by(CompoundPlate, SampleID, pert_dose, cellset, pool_id) %>% 
  dplyr::summarise(n.f1 = mean(flag1, na.rm = T),
                   n.f2 = mean(flag2, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(pmax(n.f1, n.f2) > 0.25) %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose, cellset, pool_id) %>% 
  dplyr::mutate(outlier = TRUE)



LFC.collapsed %<>% 
  dplyr::left_join(outlier.trt.pools) %>% 
  dplyr::mutate(outlier = !is.na(outlier))


LFC %<>%  
  dplyr::left_join(inst_meta) %>%
  dplyr::left_join(outlier.trt.pools) %>% 
  dplyr::distinct(prism_replicate, pert_well, analyte_id, cellset, screen, LFC, LFC_corrected, PASS, outlier) %>% 
  dplyr::mutate(outlier = !is.na(outlier))


outlier.trt.pools %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::mutate(outlier_type = "trt") %>% 
  dplyr::bind_rows(outlier.nc.pools %>%
                     dplyr::left_join(inst_meta) %>% 
                     dplyr::mutate(outlier_type = "nc")) %>% 
  dplyr::distinct(screen, CompoundPlate, prism_replicate, pert_well, pool_id, outlier_type) %>% 
  write_csv("data/25Q2/PRISMOncologyReferenceFilteredPools.csv")


# -------
# FITTING DOSE RESPONSE CURVES 
# ------

LFC.FILTERED <- LFC %>%
  dplyr::left_join(inst_meta) %>%
  dplyr::left_join(analyte_meta) %>% 
  dplyr::filter(is.finite(LFC_corrected), PASS, !outlier, pool_id != "CTLBC", pert_type == "trt_cp") %>%
  dplyr::group_by(CompoundPlate, SampleID, pert_dose, analyte_id, depmap_id, pool_id, cellset) %>%
  dplyr::filter(abs(pmin(2^LFC_corrected, 1.5) - median(pmin(2^LFC_corrected, 1.5))) < 0.75) %>%
  dplyr::group_by(CompoundPlate, SampleID, analyte_id, cellset) %>%
  dplyr::filter(length(unique(pert_dose)) > 4, n() > 12) %>%
  dplyr::ungroup()

CompoundPlates <- dplyr::distinct(LFC.FILTERED, SampleID, CompoundPlate) %>% 
  dplyr::count(CompoundPlate) %>% 
  dplyr::arrange(n)

f <- function(df){
  df %>% 
    dplyr::group_by(SampleID, CompoundPlate, analyte_id, depmap_id, pool_id, cellset, screen) %>%    
    dplyr::summarise(get_best_fit(pmin(2^LFC_corrected,1.5), pert_dose)) %>%
    dplyr::ungroup()
}

DRC <- list()
time = Sys.time()
for(ix in 1:nrow(CompoundPlates)){
  
  print(CompoundPlates[ix, ])
  print(Sys.time() - time)
  time = Sys.time()
  
  DRC[[ix]] <- CompoundPlates[ix, ] %>% 
    dplyr::select(-n) %>% 
    dplyr::left_join(LFC.FILTERED) %>% 
    dplyr::group_split(SampleID) %>% 
    parallel::mclapply(f, mc.cores = parallel::detectCores() - 1) %>%
    dplyr::bind_rows()
}

DRC <- bind_rows(DRC)

DRC %<>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::mutate(note = ifelse(is.na(note), "", note)) %>% 
  dplyr::filter(note != "REMOVE") %>% 
  dplyr::group_by(CompoundPlate, SampleID, depmap_id) %>% 
  dplyr::arrange(!is.finite(auc), note == "DEPRIORITIZE", desc(cellset)) %>% 
  dplyr::mutate(priority = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-note)


DRC %>% 
  write_csv("data/25Q2/PRISMOncologyReferenceDoseResponseParameters.csv")

# ----
# CALCULATE THE FITTED LFC VALUES
# ----

LFC.fitted <- inst_meta %>% 
  dplyr::filter(pert_type == "trt_cp") %>% 
  dplyr::distinct(screen, CompoundPlate, SampleID, pert_dose, pert_dose_unit) %>%  
  dplyr::inner_join(DRC) %>% 
  dplyr::filter(successful_fit) %>% 
  dplyr::mutate(inflection = as.numeric(inflection)) %>% 
  # dplyr::left_join(analyte_meta) %>% 
  dplyr::mutate(LFC_fitted = log2(lower_limit + (upper_limit - lower_limit) / (1 + (pert_dose / inflection)^-slope))) %>% 
  dplyr::distinct(screen, CompoundPlate, SampleID, pert_dose, pert_dose_unit, cellset, depmap_id, pool_id, LFC_fitted, priority)

LFC.collapsed %<>% 
  dplyr::full_join(LFC.fitted) 

LFC.collapsed %>%  
  dplyr::distinct(screen, CompoundPlate, SampleID, pert_dose, pert_dose_unit, cellset, pool_id, depmap_id, LFC, LFC_fitted, outlier, priority) %>%
  write_csv("data/25Q2/PRISMOncologyReferenceLFCCollapsed.csv")

# -----
# WRITE PORTAL FILES
# -----

DRC <- data.table::fread("data/25Q2/PRISMOncologyReferenceDoseResponseParameters.csv")

DRC %>% 
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundPlate, Prioritized)) %>% 
  dplyr::filter(priority == 1, Prioritized, successful_fit) %>% 
  dplyr::rename(ModelID = depmap_id, EC50 = inflection, LowerAsymptote = lower_limit, UpperAsymptote = upper_limit, Slope = slope) %>% 
  dplyr::distinct(ModelID, SampleID, CompoundPlate, EC50, LowerAsymptote, UpperAsymptote, Slope) %>% 
  write_csv("data/25Q2/PRISMOncologyReferenceResponseCurves.csv")

DRC %>% 
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundPlate, Prioritized)) %>% 
  dplyr::filter(priority == 1, Prioritized, successful_fit) %>% 
  reshape2::acast(depmap_id ~ SampleID, value.var = "log2_auc") %>% 
  write.csv("data/25Q2/PRISMOncologyReferenceLog2AUCMatrix.csv")


LFC.collapsed <- data.table::fread("data/25Q2/PRISMOncologyReferenceLFCCollapsed.csv")

LFC.collapsed %>%
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate, Prioritized)) %>% 
  dplyr::filter(priority == 1, Prioritized, !outlier) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, LFC_fitted, CompoundName) %>% 
  dplyr::mutate(Label = paste0(CompoundName, " (", SampleID, ") @", pert_dose, " ", pert_dose_unit)) %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit, depmap_id, LFC_fitted) %>% 
  reshape2::acast(depmap_id ~ Label, value.var = "LFC_fitted") %>% 
  write.csv("data/25Q2/PRISMOncologyReferenceLog2ViabilityCollapsedMatrix.csv")


LFC.collapsed %>%
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate, Prioritized)) %>% 
  dplyr::filter(priority == 1, Prioritized, !outlier) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, LFC_fitted, CompoundName) %>% 
  dplyr::mutate(Label = paste0(CompoundName, " (", SampleID, ") @", pert_dose, " ", pert_dose_unit)) %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit) %>% 
  write_csv("data/25Q2/PRISMOncologyReferenceLog2ViabilityCollapsedConditions.csv")



LFC <- data.table::fread("data/25Q2/PRISMOncologyReferenceLFC.csv")


LFC <- LFC.collapsed %>%
  dplyr::left_join(compound_annotations %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate, Prioritized)) %>% 
  dplyr::filter(priority == 1, Prioritized, !outlier) %>% 
  dplyr::distinct(screen, CompoundPlate, SampleID, pert_dose, pert_dose_unit,
                  depmap_id, pool_id, cellset) %>%
  dplyr::left_join(inst_meta) %>%
  dplyr::left_join(analyte_meta) %>% 
  dplyr::distinct(prism_replicate, pert_well, analyte_id, cellset, screen, 
                  depmap_id, SampleID, pert_dose, pert_dose_unit, CompoundPlate) %>% 
  dplyr::left_join(LFC) %>% 
  dplyr::filter(is.finite(LFC_corrected), PASS) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, CompoundPlate, LFC_corrected) %>% 
  dplyr::group_by(depmap_id, SampleID, pert_dose, pert_dose_unit, CompoundPlate) %>% 
  dplyr::mutate(Replicate = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) 

LFC.conditions <- LFC %>% 
  dplyr::distinct(SampleID, Dose, DoseUnit, CompoundPlate, Replicate) %>% 
  dplyr::arrange(SampleID, CompoundPlate, Dose, DoseUnit, Replicate) %>% 
  dplyr::mutate(Label = 1:n() - 1) 


LFC %>% 
  dplyr::left_join(LFC.conditions) %>% 
  dplyr::mutate(viability = pmin(2^LFC_corrected, 1.5)) %>% 
  reshape2::acast(depmap_id ~ Label, value.var = "viability") %>%
  write.csv("data/25Q2/PRISMOncologyReferenceViabilityMatrix.csv")

LFC.conditions %>% 
  write_csv("data/25Q2/PRISMOncologyReferenceViabilityConditions.csv")


conf.pools <- analyte_meta %>% 
  dplyr::filter(pool_id != "CTLBC", !is.na(depmap_id), !is.na(pool_id), is.na(note)) %>% 
  dplyr::mutate(ix = 1) %>%
  reshape2::acast(depmap_id ~ pool_id + cellset + screen, value.var = "ix", fill = 0) 

conf.pools <- conf.pools[, !duplicated(t(conf.pools))]

QC <- data.table::fread("data/25Q2/PRISMOncologyReferenceQCTable.csv")

conf.QC <- QC %>% 
  dplyr::inner_join(LFC.collapsed %>% 
                      dplyr::filter(priority == 1) %>% 
                      dplyr::distinct(screen, cellset, pool_id, depmap_id)) %>% 
  dplyr::distinct(depmap_id, NC.median, NC.mad, PC.median, PC.mad, DR, SSMD) %>% 
  tidyr::pivot_longer(cols = c(2:7)) %>% 
  dplyr::filter(is.finite(value)) %>% 
  reshape2::acast(depmap_id ~ name, fun.aggregate = median)

conf.pools %>% 
  reshape2::melt() %>% 
  dplyr::bind_rows(conf.QC %>% 
                     reshape2::melt()) %>% 
  reshape2::acast(Var1 ~ Var2)  %>% 
  write.csv("data/25Q2/PRISMOncologyReferenceConfounderMatrix.csv")