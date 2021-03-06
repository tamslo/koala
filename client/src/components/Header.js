import React from "react";
import styled from "styled-components";
import AppBar from "@material-ui/core/AppBar";
import Toolbar from "@material-ui/core/Toolbar";
import Typography from "@material-ui/core/Typography";
import logo from "../assets/logo.png";

export default props => {
  return (
    <AppBar position="static" color="primary">
      <Toolbar>
        <Image src={logo} alt="Koala logo" />
        <Typography variant="title" color="inherit">
          Koala RNA-Seq Platform
        </Typography>
      </Toolbar>
    </AppBar>
  );
};

const Image = styled.img`
  width: 60px;
  height: 60px;
  padding-top: 12px;
  padding-bottom: 12px;
  margin-right: 16px;
`;
