import * as types from "../actions/actionTypes";
import datasets from "./datasets";
import experiments from "./experiments";

const initialState = { context: null };

export default (state = initialState, action = {}) => {
  const { context } = state;
  if (action.type.includes("DATASET")) {
    return {
      context: { ...context, datasets: datasets(context.datasets, action) }
    };
  } else if (action.type.includes("EXPERIMENT")) {
    return {
      context: {
        ...context,
        experiments: experiments(context.experiments, action)
      }
    };
  } else {
    switch (action.type) {
      case types.FETCH_CONTEXT:
        return { context: action.context };
      default:
        return state;
    }
  }
};
