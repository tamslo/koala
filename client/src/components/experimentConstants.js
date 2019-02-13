export const ALIGNMENT = "alignment";
export const ALIGNMENT_FILTERING = "alignment_filtering";
export const VARIANT_CALLING = "variant_calling";
export const VARIANT_FILTERING = "variant_filtering";

export const serviceTypes = {
  [ALIGNMENT]: "aligner",
  [ALIGNMENT_FILTERING]: "alignment_filter",
  [VARIANT_CALLING]: "variant_caller",
  [VARIANT_FILTERING]: "variant_filter"
};

export const displayNames = {
  dataset: "Data set",
  reference: "Reference genome"
};
