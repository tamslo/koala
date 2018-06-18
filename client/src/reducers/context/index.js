import * as types from "../../actions/actionTypes";
import { addDataset } from "./datasets";
import {
  updateExperiment,
  updateExperiments,
  deleteFromExperiments
} from "./experiments";

const initialState = null;

export default (state = initialState, action = {}) => {
  switch (action.type) {
    case types.FETCH_CONTEXT:
      return action.context;
    case types.ADDING_DATASET:
      return { ...state, datasetLoading: true };
    case types.ADDED_DATASET:
      return {
        ...state,
        datasetLoading: false,
        datasets: addDataset(state.datasets, action.dataset)
      };
    case types.ADD_EXPERIMENT:
      return updateExperiments(state, action.experiment);
    case types.UPDATE_EXPERIMENT:
      const experiment = updateExperiment(state, action.experiment);
      return updateExperiments(state, experiment);
    case types.DELETE_EXPERIMENT:
      return deleteFromExperiments(state, action.experiment);
    default:
      return state;
  }
};
