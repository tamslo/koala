import * as types from "../../actions/actionTypes";

const initialState = {};

export default (state = initialState, action = {}) => {
  switch (action.type) {
    case types.FETCH_CONTEXT:
      return action.context.datasets;
    case types.ADDING_DATASET:
      return { ...state, areLoading: true };
    case types.ADDED_DATASET:
      const { dataset } = action;
      return {
        ...state,
        areLoading: false,
        [dataset.id]: dataset
      };
    default:
      return state;
  }
};
