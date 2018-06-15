const hostname = window && window.location && window.location.hostname;
export const SERVER_URL = "http://" + hostname + ":5000";

export const getRequest = (route, json = true) => {
  const request = new Request(`${SERVER_URL}${route}`);
  return handleFetch(request, json);
};

export const postRequest = (route, body) => {
  return bodyRequest(route, body, "POST");
};

export const putRequest = (route, body) => {
  return bodyRequest(route, body, "PUT");
};

export const deleteRequest = route => {
  const request = new Request(`${SERVER_URL}${route}`, { method: "DELETE" });
  return handleFetch(request);
};

const bodyRequest = (route, body, method) => {
  const request = new Request(`${SERVER_URL}${route}`, {
    headers: {
      "Content-Type": "application/json",
      "Access-Control-Allow-Origin": "*"
    },
    method,
    body: JSON.stringify(body),
    mode: "cors",
    timeout: 0
  });

  return handleFetch(request);
};

const handleFetch = (request, json = true) => {
  return fetch(request)
    .then(response => {
      if (!response.ok) {
        throw new Error(`Server returned ${response.status}`);
      }
      return json ? response.json() : response;
    })
    .then(json => json)
    .catch(error => {
      return { error };
    });
};
