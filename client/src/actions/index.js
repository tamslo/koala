import * as types from "./actionTypes";
import { getRequest, postRequest, deleteRequest } from "./_request";

export const fetchContext = () => {
  return dispatch => {
    getRequest("/context").then(context => {
      dispatch({
        type: types.FETCH_CONTEXT,
        context
      });
    });
  };
};

export const deleteExperiment = id => {
  return dispatch => {
    deleteRequest("/experiment/" + id).then(experiment => {
      dispatch({
        type: types.DELETE_EXPERIMENT,
        experiment
      });
    });
  };
};

export const runExperiment = params => {
  return dispatch => {
    postRequest("/experiment", params).then(experiment => {
      dispatch({
        type: types.ADD_EXPERIMENT,
        experiment
      });
      handleExperimentRun(experiment, dispatch);
    });
  };
};

const handleExperimentRun = (experiment, dispatch) => {
  getRequest("/data/" + experiment.id)
    .then(experiment => {
      return updateExperiment(experiment, dispatch);
    })
    .then(experiment => {
      // TODO alignment
      return experiment;
    })
    .then(experiment => {
      return getRequest("/done/" + experiment.id);
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
