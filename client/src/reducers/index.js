import * as types from "../actions/actionTypes";
import datasets from "./datasets";
import experiments from "./experiments";

const initialState = { context: null, serverUnresponsive: false };

export default (state = initialState, action = {}) => {
  const { context } = state;
  if (action.type.includes("DATASET")) {
    return {
      ...state,
      context: { ...context, datasets: datasets(context.datasets, action) }
    };
  } else if (action.type.includes("EXPERIMENT")) {
    return {
      ...state,
      context: {
        ...context,
        experiments: experiments(context.experiments, action)
      }
    };
  } else {
    switch (action.type) {
      case types.FETCH_CONTEXT:
        return { ...state, context: action.context };
      case types.SERVER_UNRESPONSIVE:
        return { ...state, serverUnresponsive: true };
      default:
        return state;
    }
  }
};
