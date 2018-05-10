import React from "react";
import Header from "./components/Header";
import Content from "./components/Content";

export default props => {
  return (
    <div className="app">
      <Header />
      <Content />
    </div>
  );
};
