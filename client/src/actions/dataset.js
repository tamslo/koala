import * as types from "./actionTypes";
import { postRequest } from "../api";

export const addDataset = dataset => {
  return dispatch => {
    dispatch({ type: types.ADDING_DATASET });

    let data = new FormData();
    Object.keys(dataset.data).forEach(fileKey =>
      data.append(fileKey, dataset.data[fileKey])
    );
    data.append("json", JSON.stringify(dataset));

    postRequest("/dataset", data, false)
      .then(dataset => {
        if (!dataset.error) {
          dispatch({
            type: types.ADDED_DATASET,
            dataset
          });
        } else {
          dispatch({ type: types.DATASET_ERROR });
        }
      })
      .catch(error => console.error(error));
  };
};
