import * as types from "./actionTypes";
import { getJson, postRequest, deleteRequest, putRequest } from "../api";
import constants from "../constants.json";

export const deleteExperiment = id => {
  return dispatch => {
    deleteRequest("/experiment?id=" + id)
      .then(experiment => {
        dispatch({
          type: types.DELETE_EXPERIMENT,
          experiment
        });
      })
      .catch(error => console.log(error));
  };
};

export const addExperiment = params => {
  return dispatch => {
    postRequest("/experiment", params)
      .then(experiment => {
        if (!experiment.error) {
          dispatch({
            type: types.ADD_EXPERIMENT,
            experiment
          });
        }
      })
      .catch(error => ({ error }));
  };
};

export const retryExperiment = experiment => {
  return dispatch => {
    const params = { ...experiment, interrupted: false };
    putRequest("/experiment", params).then(experiment => {
      dispatch({
        type: types.ADD_EXPERIMENT,
        experiment
      });
    });
  };
};

export const runExperiment = id => {
  return dispatch => {
    dispatch({
      type: types.RUN_EXPERIMENT,
      id
    });

    // Chain execution steps
    let latestExperiment = null;
    const { actions } = constants;

    getJson(`/execute?action=${actions.DATASET}&experiment=${id}`)
      .then(experiment => {
        latestExperiment = experiment;
        updateExperiment(experiment, dispatch);
        return getJson(`/execute?action=${actions.ALIGNMENT}&experiment=${id}`);
      })
      .then(experiment => {
        latestExperiment = experiment;
        updateExperiment(experiment, dispatch);
        return getJson("/done?experiment=" + id);
      })
      .then(experiment => {
        latestExperiment = experiment;
        updateExperiment(experiment, dispatch);
        dispatch({
          type: types.EXPERIMENT_DONE
        });
      })
      .catch(error => {
        const experiment = { ...latestExperiment, id, error };
        dispatch({
          type: types.UPDATE_EXPERIMENT,
          experiment
        });
        dispatch({
          type: types.EXPERIMENT_DONE
        });
      });
  };
};

const updateExperiment = (experiment, dispatch) => {
  if (experiment.error) {
    throw new Error(experiment.error);
  }

  dispatch({
    type: types.UPDATE_EXPERIMENT,
    experiment
  });

  return experiment;
};
