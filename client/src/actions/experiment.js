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

export const runExperiment = id => {
  return dispatch => {
    dispatch({
      type: types.RUN_EXPERIMENT,
      id
    });

    let result = {};
    getJson(`/execute?experiment=${id}`)
      .then(experiment => {
        result = experiment;
        if (experiment.error) {
          throw new Error(experiment.error);
        }
      })
      .catch(error => {
        const experiment = { ...result, id, error };
        result = experiment;
      })
      .finally(() => {
        dispatch({
          type: types.EXPERIMENT_DONE,
          experiment: result
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
