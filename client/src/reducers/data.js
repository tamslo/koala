import * as types from "../actions/actionTypes";

const initialState = {
  isLoading: false,
  content: null
};

export default (state = initialState, action = {}) => {
  switch (action.type) {
    case types.FETCH_DATA:
      return { ...state, isLoading: true };
    default:
      return state;
  }
};
