import * as types from "./actionTypes";
import { getJson, postRequest, deleteRequest, putRequest } from "../api";
import constants from "../constants";

export const deleteExperiment = id => {
  return dispatch => {
    deleteRequest("/experiment?id=" + id)
      .then(experiment => {
        dispatch({
          type: types.DELETE_EXPERIMENT,
          experiment
        });
      })
      .catch(error => console.error(error));
  };
};

export const addExperiment = params => {
  return dispatch => {
    postRequest("/experiment", params)
      .then(experiments => {
        dispatch({
          type: types.ADD_EXPERIMENTS,
          experiments
        });
      })
      .catch(error => ({ error }));
  };
};

export const updateRunningExperiment = () => {
  return dispatch => {
    getJson("/running")
      .then(
        experiment =>
          experiment.id &&
          dispatch({
            type: types.UPDATE_EXPERIMENT,
            experiment
          })
      )
      .catch(error => {
        console.error(error);
        if (error.toString() === "TypeError: Failed to fetch") {
          dispatch({
            type: types.SERVER_UNRESPONSIVE
          });
        }
      });
  };
};

export const retryExperiment = experiment => {
  return dispatch => {
    const params = resetStatus(experiment);
    putRequest("/experiment", params)
      .then(experiment => {
        dispatch({
          type: types.OVERWRITE_EXPERIMENT,
          experiment
        });
      })
      .catch(error => console.error(error));
  };
};

const resetStatus = experiment => {
  const pipeline = Object.keys(experiment.pipeline).reduce(
    (resettedPipeline, action) => {
      return {
        ...resettedPipeline,
        [action]: { id: experiment.pipeline[action].id }
      };
    },
    {}
  );
  return { ...experiment, status: constants.experiment.WAITING, pipeline };
};
