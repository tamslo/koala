import React, { Component } from "react";
import Header from "./components/Header";
import Content from "./components/Content";

export default class extends Component {
  render() {
    return (
      <div className="app">
        <Header />
        <Content />
      </div>
    );
  }
}
