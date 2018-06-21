import React from "react";
import withStore from "./hocs/withStore";
import withJobRunner from "./hocs/withJobRunner";
import Header from "./components/Header";
import Content from "./components/Content";

const Body = withStore(withJobRunner(Content));

export default props => {
  return (
    <div className="app">
      <Header />
      <Body />
    </div>
  );
};
