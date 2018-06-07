import * as types from "./actionTypes";
import { postRequest, deleteRequest } from "../request";

export const addDataset = dataset => {
  return dispatch => {
    postRequest("/dataset", dataset).then(dataset => {
      dispatch({
        type: types.ADD_DATASET,
        dataset
      });
    });
  };
};

export const deleteDataset = id => {
  return dispatch => {
    deleteRequest("/dataset?id=" + id).then(dataset => {
      dispatch({
        type: types.DELETE_DATASET,
        dataset
      });
    });
  };
};
