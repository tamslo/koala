import React, { Component } from "react";

let refreshInterval;

export default WrappedComponent => {
  return class extends Component {
    componentDidUpdate() {
      const { jobs, runExperiment } = this.props;
      if (jobs.running === null && jobs.waiting.length > 0) {
        runExperiment(jobs.waiting[0]);
      }
    }

    render() {
      const { jobs } = this.props;
      clearInterval(refreshInterval);

      if (jobs.running) {
        refreshInterval = setInterval(() => {
          console.log("FOO");
        }, 5000);
      }

      return <WrappedComponent {...this.props} />;
    }
  };
};
