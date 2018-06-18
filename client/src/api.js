const hostname = window && window.location && window.location.hostname;
export const SERVER_URL = "http://" + hostname + ":5000";

export const getJson = route => {
  const request = buildRequest({ route });
  return handleFetch(request);
};

export const postRequest = (route, body, jsonBody = true) => {
  const request = buildRequest({ route, method: "POST", body, jsonBody });
  return handleFetch(request);
};

export const putRequest = (route, body) => {
  const request = buildRequest({ route, method: "PUT", body });
  return handleFetch(request);
};

export const deleteRequest = route => {
  const request = buildRequest({ route, method: "DELETE" });
  return handleFetch(request);
};

const buildRequest = params => {
  const { route, method = "GET", body = null, jsonBody = true } = params;

  let headers = { "Access-Control-Allow-Origin": "*" };
  if (jsonBody) {
    headers = { ...headers, "Content-Type": "application/json" };
  }

  let options = {
    headers,
    method,
    mode: "cors",
    timeout: 0
  };

  if (body !== null) {
    options = { ...options, body: jsonBody ? JSON.stringify(body) : body };
  }

  return new Request(`${SERVER_URL}${route}`, options);
};

const handleFetch = (request, jsonResponse = true) => {
  return new Promise((resolve, reject) => {
    fetch(request)
      .then(response => {
        if (!response.ok) {
          throw new Error(`Server returned ${response.status}`);
        }
        return jsonResponse ? response.json() : response;
      })
      // in case JSON is retuned, wait for response.json() promise to resolve
      .then(json => resolve(json))
      .catch(error => reject(error));
  });
};
