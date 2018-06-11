import * as types from "./actionTypes";
import { postRequest } from "../request";

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
