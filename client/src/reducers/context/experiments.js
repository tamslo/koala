import * as types from "../../actions/actionTypes";

const initialState = {};

export default (state = initialState, action = {}) => {
  const { experiment } = action;
  switch (action.type) {
    case types.ADD_EXPERIMENT:
      return {
        ...state,
        [experiment.id]: experiment
      };
    case types.UPDATE_EXPERIMENT:
      return {
        ...state,
        [experiment.id]: { ...state[experiment.id], ...experiment }
      };
    case types.DELETE_EXPERIMENT:
      return Object.keys(state).reduce((experiments, experimentId) => {
        return experimentId === experiment.id
          ? experiments
          : { ...experiments, [experimentId]: state[experimentId] };
      }, {});
    default:
      return state;
  }
};
