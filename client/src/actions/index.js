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
      if (!experiment.isError) {
        dispatch({
          type: types.DELETE_EXPERIMENT,
          experiment
        });
      }
    });
  };
};

export const addExperiment = params => {
  return dispatch => {
    postRequest("/experiment", params).then(experiment => {
      if (!experiment.isError) {
        dispatch({
          type: types.ADD_EXPERIMENT,
          experiment
        });
      }
    });
  };
};

export const runExperiment = id => {
  return dispatch => {
    dispatch({
      type: types.RUN_EXPERIMENT,
      id
    });
    handleExperimentRun(id, dispatch)
      .then(response => {
        dispatch({
          type: types.EXPERIMENT_DONE
        });
      })
      .catch(console.error);
  };
};

const handleExperimentRun = (experimentId, dispatch) => {
  return new Promise((resolve, reject) => {
    getRequest("/data/" + experimentId)
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
        resolve(updateExperiment(experiment, dispatch));
      })
      .catch(error => reject);
  });
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
