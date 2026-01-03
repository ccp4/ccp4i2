import { useEffect, useState } from "react";
import { GeneralTable } from "../General/GeneralTable";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import $ from "jquery";

export const CCP4i2ReportTable: React.FC<CCP4i2ReportElementProps> = (
  props
) => {
  const [data, setData] = useState<any[]>([]);
  const [columns, setColumns] = useState<any[]>([]);

  useEffect(() => {
    let newData: any[] = [];
    let newColumns: any[] = [];

    $(props.item)
      .find("tr")
      .each((iRow, tableRow) => {
        var dataItem: any = { key: iRow };
        if ($(props.item).attr("transpose") === "True") {
          $(tableRow)
            .find("th")
            .each((iColumn, tableData) => {
              dataItem["col_" + iColumn] = $(tableData).text();
            });
        }
        const thCount = Object.keys(dataItem).length;
        $(tableRow)
          .find("td")
          .each((iColumn, tableData) => {
            dataItem["col_" + (thCount + iColumn - 1)] = $(tableData).text();
          });
        if ($(tableRow).find("td").length > 0) {
          newData.push(dataItem);
        }
        if ($(props.item).attr("transpose") === "False") {
          $(tableRow)
            .find("th")
            .each((iColumn, tableData) => {
              while (iColumn >= newColumns.length) {
                newColumns.push({});
              }
              newColumns[iColumn] = {
                title: (
                  <div
                    dangerouslySetInnerHTML={{ __html: tableData.innerHTML }}
                  />
                ),
                dataIndex: "col_" + iColumn,
                key: "col_" + iColumn,
                searchable: true,
                render: (text: string) => (
                  <div dangerouslySetInnerHTML={{ __html: text }} />
                ),
              };
            });
        }
      });
    if (newColumns.length === 0 && newData.length > 0) {
      for (var iCol = 0; iCol < Object.keys(newData[0]).length - 1; iCol++) {
        newColumns.push({
          title: "",
          dataIndex: "col_" + iCol,
          key: "col_" + iCol,
          searchable: true,
          render: (text: string) => (
            <div dangerouslySetInnerHTML={{ __html: text }} />
          ),
        });
      }
    }
    setData(newData);
    setColumns(newColumns);
  }, [props.job, props.item]);

  return (
    <>
      <GeneralTable
        columns={columns}
        dataSource={data}
        size="small"
        pagination={false}
        sx={{ mx: 6, my: "0.5rem", py: 0 }}
      />
    </>
  );
};
