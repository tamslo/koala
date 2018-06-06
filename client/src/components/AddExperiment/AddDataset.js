import React, { Component } from "react";
import styled from "styled-components";
import Dialog, {
  DialogActions,
  DialogTitle,
  DialogContent
} from "material-ui/Dialog";
import TextField from "material-ui/TextField";
import Button from "material-ui/Button";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = { url: "" };
  }

  handleChange = name => event => {
    this.setState({
      [name]: event.target.value
    });
  };

  render() {
    return (
      <Dialog
        open={this.props.open}
        aria-labelledby="dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="dialog-title">Add Data Set</DialogTitle>
        <StyledDialogContent>
          <DataInput
            label="Data URL"
            value={this.state.url}
            onChange={this.handleChange("url")}
            margin="normal"
          />
        </StyledDialogContent>
        <DialogActions>
          <Button onClick={this.props.cancel} color="default">
            Cancel
          </Button>
          <Button onClick={this.props.cancel} color="primary">
            Add
          </Button>
        </DialogActions>
      </Dialog>
    );
  }
}

const StyledDialogContent = styled(DialogContent)`
  display: flex;
  align-items: center;
`;

const DataInput = styled(TextField)`
-  flex-grow: 1;
-  min-width: 200px;
-  margin-right: 20px !important;
-`;
