import React, { Component } from "react";

export default WrappedComponent => {
  return class extends Component {
    componentDidUpdate() {
      const { jobs, runExperiment } = this.props;
      if (jobs.running === null && jobs.waiting.length > 0) {
        runExperiment(jobs.waiting[0]);
      }
    }

    render() {
      return <WrappedComponent {...this.props} />;
    }
  };
};
