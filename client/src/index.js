import React from "react";
import ReactDOM from "react-dom";
import { injectGlobal } from "styled-components";
import App from "./App";
import registerServiceWorker from "./registerServiceWorker";

// Theme
import Roboto from "./assets/Roboto-Regular.ttf";

// Redux
import { Provider } from "react-redux";
import { createStore, applyMiddleware, compose } from "redux";
import thunk from "redux-thunk"; // needed for async actions
import reducers from "./reducers";

const composeEnhancers = window.__REDUX_DEVTOOLS_EXTENSION_COMPOSE__ || compose;
const store = createStore(reducers, composeEnhancers(applyMiddleware(thunk)));

ReactDOM.render(
  <Provider store={store}>
    <App />
  </Provider>,
  document.getElementById("root")
);
registerServiceWorker();

injectGlobal`
  body {
    @font-face {
      font-family: "Roboto";
      font-style: normal;
      font-weight: 400;
      src: url(${Roboto});
    }

    margin: 0;
    padding: 0;
    font-family: "Roboto", sans-serif;
  }
`;
