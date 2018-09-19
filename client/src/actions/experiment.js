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

export const updateExperiment = experiment_id => {
  return dispatch => {
    getJson(`/experiment?id=${experiment_id}`).then(experiment =>
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
