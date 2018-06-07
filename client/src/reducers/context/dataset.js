export const addDataset = (datasets, dataset) => {
  return { ...datasets, [dataset.id]: dataset };
};
