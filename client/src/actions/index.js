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
    post("/experiment", experiment).then(experiment => {
      dispatch({
        type: types.ADD_EXPERIMENT,
        experiment
      });
      get(`/run?experiment=${experiment.id}`).then(report => {
        if (report.isError) {
          dispatch({
            type: types.EXPERIMENT_ERROR,
            experimentId: experiment.id,
            error: report.error.message
          });
        } else {
          dispatch({
            type: types.EXPERIMENT_DONE,
            experimentId: experiment.id,
            report
          });
        }
      });
    });
  };
};
