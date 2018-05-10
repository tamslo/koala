import * as types from "../actions/actionTypes";

const initialState = null;

export default (state = initialState, action = {}) => {
  switch (action.type) {
    case types.FETCH_CONTEXT:
      return action.context;
    case types.ADD_EXPERIMENT:
      return updateExperiments(state, { ...action.experiment, done: false });
    case types.UPDATE_EXPERIMENT:
      const experiment = updateExperiment(state, action.experiment);
      return updateExperiments(state, experiment);
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
