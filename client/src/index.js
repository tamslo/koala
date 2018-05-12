import React from "react";
import ReactDOM from "react-dom";
import { injectGlobal } from "styled-components";
import App from "./App";
import registerServiceWorker from "./registerServiceWorker";

// Theme
import { MuiThemeProvider, createMuiTheme } from "material-ui/styles";
import Roboto from "./assets/Roboto-Regular.ttf";
import { green, yellow } from "material-ui/colors";

// Redux
import { Provider } from "react-redux";
import { createStore, applyMiddleware, compose } from "redux";
import thunk from "redux-thunk"; // needed for async actions
import reducers from "./reducers";

const composeEnhancers = window.__REDUX_DEVTOOLS_EXTENSION_COMPOSE__ || compose;
const store = createStore(reducers, composeEnhancers(applyMiddleware(thunk)));

const theme = createMuiTheme({
  palette: {
    primary: { main: green[800] },
    secondary: { main: green[600] },
    warning: { main: yellow[600] }
  }
});

ReactDOM.render(
  <MuiThemeProvider theme={theme}>
    <Provider store={store}>
      <App />
    </Provider>
  </MuiThemeProvider>,
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
    background-color: #eeeeee;
  }
`;
