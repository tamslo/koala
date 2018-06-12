import * as types from "./actionTypes";
import { postRequest } from "../request";

export const addDataset = dataset => {
  return dispatch => {
    // postRequest("/dataset", dataset).then(dataset => {
    //   dispatch({
    //     type: types.ADD_DATASET,
    //     dataset
    //   });
    // });

    let data = new FormData();
    data.append("file", dataset.content);
    data.append("json", JSON.stringify(dataset));

    const request = new Request("http://localhost:5000/dataset", {
      headers: {
        "Access-Control-Allow-Origin": "*"
      },
      method: "POST",
      body: data,
      mode: "cors",
      timeout: 0
    });

    fetch(request)
      .then(response => {
        if (!response.ok) {
          throw new Error(`Server returned ${response.status}`);
        }
        return response.json();
      })
      .catch(error => {
        return { error };
      });
  };
};
