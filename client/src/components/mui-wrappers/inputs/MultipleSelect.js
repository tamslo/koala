import React from "react";
import Input from "@material-ui/core/Input";
import InputLabel from "@material-ui/core/InputLabel";
import MenuItem from "@material-ui/core/MenuItem";
import FormControl from "@material-ui/core/FormControl";
import ListItemText from "@material-ui/core/ListItemText";
import Select from "@material-ui/core/Select";
import Checkbox from "@material-ui/core/Checkbox";
import { width, marginTop, marginRight, marginBottom } from "./constants";

export default props => {
  const textFieldStyle = {
    marginTop,
    marginRight,
    marginBottom,
    width
  };
  return (
    <FormControl style={textFieldStyle}>
      <InputLabel htmlFor="select-multiple-checkbox">{props.label}</InputLabel>
      <Select
        multiple
        value={props.selected}
        onChange={props.onChange}
        input={<Input id="select-multiple-checkbox" />}
        renderValue={selected =>
          selected
            .map(item_id => {
              const item = props.items.find(item => item.id === item_id);
              return item && item.name;
            })
            .join(", ")
        }
      >
        {props.items.map(item => (
          <MenuItem key={item.id} value={item.id}>
            <Checkbox checked={itemSelected(props.selected, item.id)} />
            {item.content || <ListItemText primary={item.name} />}
          </MenuItem>
        ))}
        {props.children}
      </Select>
    </FormControl>
  );
};

const itemSelected = (selected, item_id) => {
  return typeof selected.find(id => id === item_id) !== "undefined";
};
