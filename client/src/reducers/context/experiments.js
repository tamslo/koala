export const updateExperiments = (state, experiment) => {
  return {
    ...state,
    experiments: {
      ...state.experiments,
      [experiment.id]: experiment
    }
  };
};

export const updateExperiment = (state, experiment) => {
  return { ...state.experiments[experiment.id], ...experiment };
};

export const deleteFromExperiments = (state, experiment) => {
  const experiments = Object.keys(state.experiments).reduce(
    (experiments, experimentId) => {
      return experimentId === experiment.id
        ? experiments
        : { ...experiments, [experimentId]: state.experiments[experimentId] };
    },
    {}
  );
  return { ...state, experiments };
};
