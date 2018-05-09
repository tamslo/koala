import * as types from "./actionTypes";
import { get, post } from "./_request";

export const fetchContext = () => {
  return dispatch => {
    get("/context").then(context => {
      dispatch({
        type: types.FETCH_CONTEXT,
        context
      });
    });
  };
};

export const run = experiment => {
  return dispatch => {
    dispatch({
      type: types.ADD_EXPERIMENT,
      experiment
    });
    post("/run", experiment).then(experiment => {
      dispatch({
        type: types.EXPERIMENT_DONE,
        experiment
      });
    });
  };
};
