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

export const run = params => {
  return dispatch => {
    post("/experiment", params).then(experiment => {
      dispatch({
        type: types.ADD_EXPERIMENT,
        experiment
      });
      runExperiment(experiment, dispatch);
    });
  };
};

const runExperiment = (experiment, dispatch) => {
  get(`/data?experiment=${experiment.id}`)
    .then(experiment => {
      return updateExperiment(experiment, dispatch);
    })
    .then(experiment => {
      // TODO alignment
      return experiment;
    })
    .then(experiment => {
      return get(`/done?experiment=${experiment.id}`);
    })
    .then(experiment => {
      updateExperiment(experiment, dispatch);
    })
    .catch(error => error);
};

const updateExperiment = (experiment, dispatch) => {
  dispatch({
    type: types.UPDATE_EXPERIMENT,
    experiment
  });
  if (experiment.error) {
    throw new Error(experiment);
  }
  return experiment;
};
