import List from '@mui/material/List';
import { useState, useEffect } from 'react';
import $ from 'jquery';
import { Paper } from '@mui/material';

export const CCP4i2BaseBrowser = (props) => {
    const [items, setItems] = useState([])

    useEffect(() => {
        if (props.predicate) {
            const asyncFunc = async () => {
                const paramDict = {
                    __type__: props.dbType,
                    ...props.predicate
                }
                if (props.fields) {
                    paramDict.__values__ = props.fields
                }
                const itemsResult = await fetch(`${props.urlRoot}/ModelValues?${$.param(paramDict)}`)
                    .then(response => response.json())
                setItems(props.sort ? itemsResult.results.sort(props.sort) : itemsResult.results)
                if (props.repeatLogic && props.repeatLogic(itemsResult.results)) {
                    setTimeout(async () => { await asyncFunc() }, 5000)
                }
                if (props.onLoaded) {
                    props.onLoaded(itemsResult)
                }
            }
            asyncFunc()
        }
    }, [props.predicate])

    return <>
        {props.header && <h5>{props.header}</h5>}
        <Paper>
            <div style={{ maxidth: "20%", height: "20rem", overflow: "auto" }}>
                <List dense={true}>
                    {
                        items.map(item => props.listItem(item))
                    }
                </List>
            </div>
        </Paper>
    </>
}
CCP4i2BaseBrowser.defaultProps = { predicate: {}, onSelect: () => { } }