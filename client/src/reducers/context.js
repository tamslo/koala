import * as types from "../actions/actionTypes";

const initialState = null;

export default (state = initialState, action = {}) => {
  switch (action.type) {
    case types.FETCH_CONTEXT:
      return action.context;
    case types.ADD_EXPERIMENT:
      return updateExperiments(state, action.experiment);
    case types.UPDATE_EXPERIMENT:
      const experiment = updateExperiment(state, action.experiment);
      return updateExperiments(state, experiment);
    case types.DELETE_EXPERIMENT:
      return deleteFromExperiments(state, action.experiment);
    default:
      return state;
  }
};

const updateExperiments = (state, experiment) => {
  return {
    ...state,
    experiments: {
      ...state.experiments,
      [experiment.id]: experiment
    }
  };
};

const updateExperiment = (state, experiment) => {
  return { ...state.experiments[experiment.id], ...experiment };
};

const deleteFromExperiments = (state, experiment) => {
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
