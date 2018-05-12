import * as types from "../actions/actionTypes";

const initialState = { running: null, waiting: [] };

export default (state = initialState, action = {}) => {
  switch (action.type) {
    case types.FETCH_CONTEXT:
      return {
        ...state,
        waiting: Object.keys(action.context.experiments).reduce(
          (waiting, experimentId) => {
            const experiment = action.context.experiments[experimentId];
            return experiment.done || experiment.interrupted || experiment.error
              ? waiting
              : [...waiting, experimentId];
          },
          []
        )
      };
    case types.RUN_EXPERIMENT:
      return {
        running: action.id,
        waiting: state.waiting.filter(id => id !== action.id)
      };
    case types.ADD_EXPERIMENT:
      return { ...state, waiting: [...state.waiting, action.experiment.id] };
    case types.DELETE_EXPERIMENT:
      if (state.waiting.includes(action.experiment.id)) {
        return {
          ...state,
          waiting: state.waiting.filter(id => id !== action.experiment.id)
        };
      }
      return { ...state, waiting: [...state.waiting, action.experiment.id] };
    case types.EXPERIMENT_DONE:
      return { ...state, running: null };
    default:
      return state;
  }
};
