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
          [action.experiment.id]: {
            ...action.experiment,
            done: false,
            error: false
          }
        }
      };
    case types.EXPERIMENT_DONE:
      return {
        ...state,
        experiments: {
          ...state.experiments,
          [action.report.id]: {
            ...state.experiments[action.report.id],
            done: true,
            report: action.report.content
          }
        }
      };
    case types.EXPERIMENT_ERROR:
      return {
        ...state,
        experiments: {
          ...state.experiments,
          [action.experimentId]: {
            ...state.experiments[action.experimentId],
            error: true
          }
        }
      };
    default:
      return state;
  }
};
