import React, { Component } from "react";
import { connect } from "react-redux";
import styled from "styled-components";
import {
  fetchContext,
  addExperiment,
  deleteExperiment,
  runExperiment
} from "../actions";
import Loading from "./Loading";
import Inputs from "./Inputs";
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
    const { context, addExperiment, deleteExperiment, jobs } = this.props;
    if (context === null) {
      return <Loading content={"Setting everything up..."} />;
    }

    if (context.isError) {
      return (
        <Loading
          content={context.error.message}
          error={context.isError}
          retry={this.props.fetchContext}
        />
      );
    }

    return (
      <Container className="content">
        <Inputs aligners={context.aligners} addExperiment={addExperiment} />
        <Spacer />
        {Object.keys(context.experiments).length > 0 && (
          <Experiments
            experiments={context.experiments}
            deleteExperiment={deleteExperiment}
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

const mapDispatchToProps = dispatch => {
  return {
    fetchContext: () => {
      dispatch(fetchContext());
    },
    addExperiment: params => {
      dispatch(addExperiment(params));
    },
    deleteExperiment: id => {
      dispatch(deleteExperiment(id));
    },
    runExperiment: id => {
      dispatch(runExperiment(id));
    }
  };
};

export default connect(mapStateToProps, mapDispatchToProps)(Content);
