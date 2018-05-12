const hostname = window && window.location && window.location.hostname;
const SERVER_URL = "http://" + hostname + ":5000";

export const getRequest = route => {
  const request = new Request(`${SERVER_URL}${route}`);
  return handleFetch(request);
};

export const postRequest = (route, body) => {
  const request = new Request(`${SERVER_URL}${route}`, {
    headers: {
      "Content-Type": "application/json",
      "Access-Control-Allow-Origin": "*"
    },
    method: "POST",
    body: JSON.stringify(body),
    mode: "cors",
    timeout: 0
  });

  return handleFetch(request);
};

export const deleteRequest = route => {
  const request = new Request(`${SERVER_URL}${route}`, { method: "DELETE" });
  return handleFetch(request);
};

const handleFetch = request => {
  return fetch(request)
    .then(response => {
      if (!response.ok) {
        throw new Error(`Server returned ${response.status}`);
      }
      return response.json();
    })
    .catch(error => {
      console.error(`Error in ${request.url}: ${error}`);
      return { error };
    });
};
