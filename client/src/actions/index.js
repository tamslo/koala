import * as types from "./actionTypes";
import { getJson } from "../api";

export const fetchContext = () => {
  return dispatch => {
    getJson("/context")
      .then(context => {
        dispatch({
          type: types.FETCH_CONTEXT,
          context
        });
      })
      .catch(error => console.log(error));
  };
};
