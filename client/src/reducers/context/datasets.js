import * as types from "../../actions/actionTypes";

const initialState = {};

export default (state = initialState, action = {}) => {
  switch (action.type) {
    case types.ADDING_DATASET:
      return { ...state, areLoading: true };
    case types.DATASET_ERROR:
      return { ...state, areLoading: false };
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
