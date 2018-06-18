import * as types from "./actionTypes";

export const addDataset = dataset => {
  return dispatch => {
    let data = new FormData();
    Object.keys(dataset.content).forEach(fileKey =>
      data.append(fileKey, dataset.content[fileKey])
    );
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

    dispatch({ type: types.ADDING_DATASET });

    fetch(request)
      .then(response => {
        if (!response.ok) {
          throw new Error(`Server returned ${response.status}`);
        }
        return response.json();
      })
      .then(dataset =>
        dispatch({
          type: types.ADDED_DATASET,
          dataset
        })
      )
      .catch(error => {
        return { error };
      });
  };
};
