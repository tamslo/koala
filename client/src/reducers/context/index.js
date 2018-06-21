import * as types from "../../actions/actionTypes";
import datasets from "./datasets";
import experiments from "./experiments";

const initialState = null;

export default (state = initialState, action = {}) => {
  if (action.type.includes("DATASET")) {
    return { ...state, datasets: datasets(state.datasets, action) };
  } else if (action.type.includes("EXPERIMENT")) {
    return { ...state, experiments: experiments(state.experiments, action) };
  } else {
    switch (action.type) {
      case types.FETCH_CONTEXT:
        console.log(action.context);
        return action.context;
      default:
        return state;
    }
  }
};
