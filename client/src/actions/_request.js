const hostname = window && window.location && window.location.hostname;
const SERVER_URL = "http://" + hostname + ":5000";

export const get = route => {
  return fetch(`${SERVER_URL}${route}`)
    .then(response => {
      if (!response.ok) {
        throw new Error(`Server returned ${response.status}`);
      }
      return response.json();
    })
    .then(json => json)
    .catch(error => {
      console.error(error);
      return { isError: true, error };
    });
};

export const post = (route, body) => {
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
  return fetch(request)
    .then(response => {
      if (!response.ok) {
        throw new Error(`Server returned ${response.status}`);
      }
      return response.json();
    })
    .then(json => json)
    .catch(error => {
      console.error(error);
      return { isError: true, error };
    });
};
