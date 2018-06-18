import React, { Component } from "react";
import { connect } from "react-redux";
import styled from "styled-components";
import {
  fetchContext,
  addExperiment,
  deleteExperiment,
  runExperiment,
  retryExperiment
} from "../actions";
import { addDataset } from "../actions/dataset";
import Loading from "./Loading";
import AddExperiment from "./AddExperiment";
import Experiments from "./Experiments";

class Content extends Component {
  componentWillMount() {
    this.props.fetchContext();
  }

  componentDidUpdate() {
    const { jobs, runExperiment } = this.props;
    if (jobs.running === null && jobs.waiting.length > 0) {
      runExperiment(jobs.waiting[0]);
    }
  }

  render() {
    const {
      context,
      addExperiment,
      deleteExperiment,
      retryExperiment,
      addDataset,
      jobs
    } = this.props;
    if (context === null) {
      return <Loading content={"Setting everything up..."} />;
    }

    if (context.error) {
      return (
        <Loading
          content={context.error.message}
          error={context.error}
          retry={this.props.fetchContext}
        />
      );
    }

    return (
      <Container className="content">
        <AddExperiment
          services={context.services}
          addExperiment={addExperiment}
          datasets={context.datasets}
          datasetLoading={context.datasetLoading}
          addDataset={addDataset}
        />
        <Spacer />
        {Object.keys(context.experiments).length > 0 && (
          <Experiments
            experiments={context.experiments}
            services={context.services}
            datasets={context.datasets}
            deleteExperiment={deleteExperiment}
            retryExperiment={retryExperiment}
            jobs={jobs}
          />
        )}
      </Container>
    );
  }
}

const Container = styled.div`
  padding: 32px;
`;

const Spacer = styled.div`
  height: 32px;
`;

const mapStateToProps = state => {
  return {
    context: state.context,
    jobs: state.jobs
  };
};

const actions = {
  fetchContext,
  addExperiment,
  deleteExperiment,
  runExperiment,
  retryExperiment,
  addDataset
};

export default connect(
  mapStateToProps,
  actions
)(Content);
