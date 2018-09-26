import * as types from "./actionTypes";
import { getJson, postRequest, deleteRequest, putRequest } from "../api";

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
      .then(experiment => {
        experiment.id &&
          dispatch({
            type: types.ADD_EXPERIMENT,
            experiment
          });
      })
      .catch(error => ({ error }));
  };
};

export const updateRunningExperiment = () => {
  return dispatch => {
    getJson("/running").then(
      experiment =>
        experiment.id &&
        dispatch({
          type: types.UPDATE_EXPERIMENT,
          experiment
        })
    );
  };
};

export const retryExperiment = experiment => {
  return dispatch => {
    const params = resetStatus(experiment);
    putRequest("/experiment", params).then(experiment => {
      dispatch({
        type: types.ADD_EXPERIMENT,
        experiment
      });
    });
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
  return { ...experiment, interrupted: false, error: false, pipeline };
};
