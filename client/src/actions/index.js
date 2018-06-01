import * as types from "./actionTypes";
import { getRequest, postRequest, deleteRequest, putRequest } from "../request";

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
    deleteRequest("/experiment?id=" + id).then(experiment => {
      dispatch({
        type: types.DELETE_EXPERIMENT,
        experiment
      });
    });
  };
};

export const addExperiment = params => {
  return dispatch => {
    postRequest("/experiment", params).then(experiment => {
      if (!experiment.error) {
        dispatch({
          type: types.ADD_EXPERIMENT,
          experiment
        });
      }
    });
  };
};

export const retryExperiment = experiment => {
  return dispatch => {
    const params = { ...experiment, interrupted: false };
    putRequest("/experiment", params).then(experiment => {
      if (!experiment.error) {
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

    // Chain execution steps
    let latestExperiment = null;
    getRequest("/execute?action=dataset&experiment=" + id)
      .then(experiment => {
        latestExperiment = experiment;
        updateExperiment(experiment, dispatch);
        return getRequest("/execute?action=alignment&experiment=" + id);
      })
      .then(experiment => {
        latestExperiment = experiment;
        updateExperiment(experiment, dispatch);
        return getRequest("/done?experiment=" + id);
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
