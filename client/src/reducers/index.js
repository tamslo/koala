import * as types from "../actions/actionTypes";

const initialState = null;

export default (state = initialState, action = {}) => {
  switch (action.type) {
    case types.FETCH_CONTEXT:
      return action.context;
    case types.ADD_EXPERIMENT:
      return {
        ...state,
        experiments: {
          ...state.experiments,
          [action.experiment.id]: { ...action.experiment, done: false }
        }
      };
    case types.EXPERIMENT_DONE:
      return {
        ...state,
        experiments: {
          ...state.experiments,
          [action.experiment.id]: { ...action.experiment, done: true }
        }
      };
    default:
      return state;
  }
};
