import React from "react";
import withStore from "./hocs/withStore";
import withUpdater from "./hocs/withUpdater";
import Header from "./components/Header";
import Content from "./components/Content";

const Body = withStore(withUpdater(Content));

export default props => {
  return (
    <div className="app">
      <Header />
      <Body />
    </div>
  );
};
